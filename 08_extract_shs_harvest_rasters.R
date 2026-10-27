# =============================================================================
# 08_shs_harvest_rasters.R
# Harvest rasters for Andy's ST-Sim runs: one raster per 5-yr period (or per year) for each
# Canfor harvest plan (8035 BAU, 8034 bird-constrained), on a grid covering
# the FMA + 10 km
#
# Cell value: 1 = harvested in that period (the cell's centre falls inside a
# harvested polygon), 0 = not
#
# Any raster turns block edges into cell edges, so some block shape is lost;
# how much depends on the cell size. Section 3 measures it on the real blocks
# at 25, 50 and 100 m so the cell size can be picked with numbers. 
#
# Output: step_lookup.csv (time step -> years from plan start -> calendar years)
# and an FMA mask on the same grid
#
# =============================================================================

library(sf)
library(terra)
library(dplyr)

proj_dir <- "C:/Users/hartt/Documents/Chapter 5/Analysis"
out_dir  <- file.path(proj_dir, "data/derived")
ras_dir  <- file.path(out_dir, "harvest_rasters_for_andy")
dir.create(ras_dir, showWarnings = FALSE)

# ---- Settings  --------------------------------------------------
cell_size <- 50             # m
grid_crs  <- "EPSG:3400"    # Alberta 10-TM Forest
buffer_m  <- 10000          # FMA + 10 km
year0     <- 2027           # calendar year of SHS year 0 
time_field <- "period"      # "period" = one raster per 5-yr period; "harvyear" = one raster per year

# ---- 1. Grid: FMA + 10 km, edges snapped to whole cells ---------------------------------
fma <- st_read(file.path(out_dir, "fma_boundary.gpkg"), quiet = TRUE)
fma <- st_union(st_transform(fma, grid_crs))
area_box <- ext(vect(st_buffer(fma, buffer_m)))

template <- rast(xmin = floor(xmin(area_box) / cell_size) * cell_size,
                 xmax = ceiling(xmax(area_box) / cell_size) * cell_size,
                 ymin = floor(ymin(area_box) / cell_size) * cell_size,
                 ymax = ceiling(ymax(area_box) / cell_size) * cell_size,
                 resolution = cell_size, crs = grid_crs)
cat("Grid:", nrow(template), "rows x", ncol(template), "cols at", cell_size, "m\n")

# ---- 2. Read the two harvest plans --------------------------------------
shs_bau  <- st_read(file.path(proj_dir, "data/raw/fmp8035_shs_100yr/fmp8035_shs_100yr.shp"), quiet = TRUE)
shs_bird <- st_read(file.path(proj_dir, "data/raw/fmp8034_shs_100yr/fmp8034_shs_100yr.shp"), quiet = TRUE)

shs_bau  <- st_transform(shs_bau, grid_crs)
shs_bird <- st_transform(shs_bird, grid_crs)

# Some SHS polygons are invalid (self-intersecting), which breaks st_union so repair here
cat("Invalid polygons: BAU", sum(!st_is_valid(shs_bau)), " bird", sum(!st_is_valid(shs_bird)), "\n")
shs_bau  <- st_make_valid(shs_bau)
shs_bird <- st_make_valid(shs_bird)

# Sort by the time step so output layers come out in time order
shs_bau  <- arrange(shs_bau,  .data[[time_field]])
shs_bird <- arrange(shs_bird, .data[[time_field]])

# ---- 3. Check how much block shape is lost at 25, 50 and 100 m (period 2, 8035) --------------
# 3a. Block sizes: merge touching polygons cut in the same period into blocks
p2 <- filter(shs_bau, period == 2)
p2_blocks <- st_cast(st_union(p2), "POLYGON")
block_ha  <- as.numeric(st_area(p2_blocks)) / 10000
cat("\nPeriod 2 (BAU):", length(block_ha), "blocks; size quantiles (ha):\n")
print(round(quantile(block_ha, c(0.05, 0.25, 0.5, 0.75, 0.95)), 2))
cat("Blocks < 1 ha:", sum(block_ha < 1), " (", round(100 * sum(block_ha[block_ha < 1]) / sum(block_ha), 1), "% of area )\n")
cat("Blocks < 0.25 ha (smaller than one 50 m cell):", sum(block_ha < 0.25), "\n")

# 3b. Shape test on a 4 x 4 km tile around the biggest period-2 block.
# "Truth" = blocks at 5 m. Each coarser grid is compared with it, cell by cell.
centre <- st_coordinates(st_centroid(p2_blocks[which.max(block_ha)]))
tx <- floor(centre[1] / 100) * 100              # snap to 100 m so all grids line up
ty <- floor(centre[2] / 100) * 100
tile <- ext(tx - 2000, tx + 2000, ty - 2000, ty + 2000)
p2_tile <- vect(st_crop(p2, st_bbox(c(xmin = tx - 2000, xmax = tx + 2000, ymin = ty - 2000, ymax = ty + 2000),
                                    crs = st_crs(p2))))

truth5 <- rasterize(p2_tile, rast(tile, resolution = 5,   crs = grid_crs), background = 0)
r25    <- rasterize(p2_tile, rast(tile, resolution = 25,  crs = grid_crs), background = 0)
r50    <- rasterize(p2_tile, rast(tile, resolution = 50,  crs = grid_crs), background = 0)
r100   <- rasterize(p2_tile, rast(tile, resolution = 100, crs = grid_crs), background = 0)

# Expand each coarse grid back up to 5 m so it lines up with the truth
b25  <- disagg(r25,  fact = 5)
b50  <- disagg(r50,  fact = 10)
b100 <- disagg(r100, fact = 20)

true_ha <- global(truth5, "sum")[1, 1] * 25 / 10000
shape_test <- data.frame(
  cell_m        = c(25, 50, 100),
  area_ha       = round(c(global(b25, "sum")[1, 1], global(b50, "sum")[1, 1],
                          global(b100, "sum")[1, 1]) * 25 / 10000, 1),
  pct_misplaced = round(100 * c(global(abs(b25  - truth5), "sum")[1, 1],
                                global(abs(b50  - truth5), "sum")[1, 1],
                                global(abs(b100 - truth5), "sum")[1, 1]) * 25 / 10000 / true_ha, 1))
cat("\nShape test, 4 x 4 km tile, period 2. True harvested area:", round(true_ha, 1), "ha\n")
print(shape_test)
# pct_misplaced = area wrongly marked cut + area wrongly marked uncut, as % of the
# true harvested area. Pick the coarsest cell size where this is small enough

png(file.path(proj_dir, "outputs", "harvest_cell_size_test.png"), width = 2000, height = 700, res = 150)
par(mfrow = c(1, 3))
plot(r25,  col = c("grey95", "darkorange"), legend = FALSE, main = "25 m cells (outline = true blocks)"); plot(p2_tile, add = TRUE, border = "black", lwd = 0.6)
plot(r50,  col = c("grey95", "darkorange"), legend = FALSE, main = "50 m cells");  plot(p2_tile, add = TRUE, border = "black", lwd = 0.6)
plot(r100, col = c("grey95", "darkorange"), legend = FALSE, main = "100 m cells"); plot(p2_tile, add = TRUE, border = "black", lwd = 0.6)
dev.off()

# ---- 4. Time-step lookup table --------------------------------------------------------------
steps_bau <- shs_bau %>% st_drop_geometry() %>% group_by(step = .data[[time_field]]) %>%
  summarise(first_year = min(harvyear), last_year = max(harvyear), ha = sum(areaha))
steps_bird <- shs_bird %>% st_drop_geometry() %>% group_by(step = .data[[time_field]]) %>%
  summarise(first_year = min(harvyear), last_year = max(harvyear), ha = sum(areaha))

steps_bau$layer_name  <- paste0(time_field, sprintf("%03d", steps_bau$step),  "_", year0 + steps_bau$first_year,  "-", year0 + steps_bau$last_year)
steps_bird$layer_name <- paste0(time_field, sprintf("%03d", steps_bird$step), "_", year0 + steps_bird$first_year, "-", year0 + steps_bird$last_year)
print(steps_bau)

# ---- 5. Rasterize and write: one 0/1 raster per time step ---------------------------------
# One layer at a time, written straight to disk - holding all layers in memory at once
# needs ~3 GB for 20 periods and ~14 GB for 100 years at 50 m so use loop 
dir.create(file.path(ras_dir, "8035_BAU"),  showWarnings = FALSE)
dir.create(file.path(ras_dir, "8034_bird"), showWarnings = FALSE)
steps_bau$ras_ha  <- NA
steps_bird$ras_ha <- NA

for (i in 1:nrow(steps_bau)) {
  this_step <- shs_bau[shs_bau[[time_field]] == steps_bau$step[i], ]
  r <- rasterize(vect(this_step), template, background = 0,
                 filename = file.path(ras_dir, "8035_BAU", paste0("harvest_8035_", steps_bau$layer_name[i], ".tif")),
                 datatype = "INT1U", overwrite = TRUE)
  steps_bau$ras_ha[i] <- global(r, "sum")[1, 1] * cell_size^2 / 10000
}

for (i in 1:nrow(steps_bird)) {
  this_step <- shs_bird[shs_bird[[time_field]] == steps_bird$step[i], ]
  r <- rasterize(vect(this_step), template, background = 0,
                 filename = file.path(ras_dir, "8034_bird", paste0("harvest_8034_", steps_bird$layer_name[i], ".tif")),
                 datatype = "INT1U", overwrite = TRUE)
  steps_bird$ras_ha[i] <- global(r, "sum")[1, 1] * cell_size^2 / 10000
}

# ---- 6. Check: harvested area per time step, rasters vs shapefile ---------------------------
print(data.frame(step = steps_bau$step,  shp_ha = round(steps_bau$ha),  ras_ha = round(steps_bau$ras_ha)))
print(data.frame(step = steps_bird$step, shp_ha = round(steps_bird$ha), ras_ha = round(steps_bird$ras_ha)))
# ras_ha should be close to shp_ha
write.csv(steps_bau,  file.path(ras_dir, "step_lookup_8035.csv"), row.names = FALSE)
write.csv(steps_bird, file.path(ras_dir, "step_lookup_8034.csv"), row.names = FALSE)

# ---- 7. Extras -------------------------------------------------------------------------------
# FMA mask on the same grid (1 = FMA), so Andy can report inside the FMA
fma_mask <- rasterize(vect(fma), template, background = 0)
writeRaster(fma_mask, file.path(ras_dir, "fma_mask.tif"), overwrite = TRUE, datatype = "INT1U")

# Quick look: first time step of each plan
first_bau  <- rast(file.path(ras_dir, "8035_BAU",  paste0("harvest_8035_", steps_bau$layer_name[1],  ".tif")))
first_bird <- rast(file.path(ras_dir, "8034_bird", paste0("harvest_8034_", steps_bird$layer_name[1], ".tif")))
png(file.path(proj_dir, "outputs", "harvest_rasters_first_step.png"), width = 1800, height = 1000, res = 150)
par(mfrow = c(1, 2))
plot(first_bau,  main = paste("8035 BAU,",  steps_bau$layer_name[1]))
plot(fma, add = TRUE, border = "grey40")
plot(first_bird, main = paste("8034 bird,", steps_bird$layer_name[1]))
plot(fma, add = TRUE, border = "grey40")
dev.off()

cat("\nFiles written to:", ras_dir, "\n")