# =========================
# Trädkron-segmentering – DUAL PASS (global + fine)
# terra mot disk + robust CRS (EPSG:3006) + mcws via raster/sp
# =========================

# Rensa minnet (valfritt)
rm(list = ls())
gc()

# Paket
library(terra)        # raster- och vektordata
library(sf)           # geometrier + CRS
library(ForestTools)  # vwf + mcws
library(dplyr)        # datahantering
library(igraph)       # glesning via komponenter
library(units)        # enheter (m^2)
library(raster)       # RasterLayer + sp-konvertering (för mcws)
library(sp)

# ====== INDATA ======
chm_path <- "C:/Path/To/Your/CHM.tif"

# ====== terra: jobba mot disk ======
.cache_root <- "C:/.Lunds Universitet/MasterThesis"   # ändra om du vill
dir.create(.cache_root, showWarnings = FALSE, recursive = TRUE)
.cache_dir <- file.path(.cache_root, "_terra_cache")
dir.create(.cache_dir, showWarnings = FALSE, recursive = TRUE)

terraOptions(
  todisk   = TRUE,          # streama blockvis till disk
  memfrac  = 0.6,           # lämna RAM-marginal
  progress = 1,             # progress i konsolen
  tempdir  = .cache_dir     # terra-tempfiler
)

# komprimerad/BigTIFF-output för mellansteg
.wopt <- list(gdal = "COMPRESS=LZW,BIGTIFF=YES")

# ====== 1) Läs CHM och skala till meter vid behov ======
chm_raw <- rast(chm_path)
mx <- global(chm_raw, "max", na.rm = TRUE)[1]
chm0 <- if (mx > 60) chm_raw / 10 else chm_raw

# --- Projektera CHM till EPSG:3006 (meter) på disk ---
target_epsg <- "EPSG:3006"
chm <- terra::project(
  chm0, target_epsg,
  filename  = file.path(.cache_dir, "chm_epsg3006.tif"),
  overwrite = TRUE, wopt = .wopt
)

# ====== 2) Förbehandling till fil ======
# Svag 3x3-glättning
k <- matrix(1, 3, 3) / 9
chm_s <- focal(
  chm, w = k, fun = mean, na.policy = "omit", pad = TRUE,
  filename  = file.path(.cache_dir, "chm_s.tif"),
  overwrite = TRUE, wopt = .wopt
)

# Maskera <2 m
chm_s <- mask(
  chm_s, chm_s > 2,
  filename  = file.path(.cache_dir, "chm_s_masked.tif"),
  overwrite = TRUE, wopt = .wopt
)

# ====== 3) Treetops – GLOBAL + FINE ======
# 3a) global fönster (lite större)
tops_global <- ForestTools::vwf(
  chm_s,
  winFun    = function(h_m) 0.045 * h_m + 0.9,
  minHeight = 3
)
tops_g <- st_as_sf(tops_global)
if (is.na(st_crs(tops_g))) st_crs(tops_g) <- st_crs(chm_s)
tops_g$source <- "global"

# 3b) fine fönster (mindre – bättre separation i tätt)
tops_fine <- ForestTools::vwf(
  chm_s,
  winFun    = function(h_m) 0.035 * h_m + 0.8,
  minHeight = 3
)
tops_f <- st_as_sf(tops_fine)
if (is.na(st_crs(tops_f))) st_crs(tops_f) <- st_crs(chm_s)
if (st_crs(tops_f) != st_crs(tops_g)) tops_f <- st_transform(tops_f, st_crs(tops_g))
tops_f$source <- "fine"

# Slå ihop toppar
toppar_sf <- rbind(tops_g, tops_f)

# ====== 3c) Glesning 0.8 m med prio: fine > global, därefter CHM-höjd ======
# Höjd vid topp (sekundär prioritet)
toppar_sf$h_chm <- terra::extract(chm_s, terra::vect(toppar_sf))[, 2]  # ev. varning är OK

# Kluster via buffert (meter, eftersom vi nu är i EPSG:3006)
buf  <- sf::st_buffer(toppar_sf, 0.8)
grp  <- sf::st_intersects(buf)
memb <- igraph::components(igraph::graph_from_adj_list(grp))$membership
toppar_sf$grp  <- memb
toppar_sf$prio <- ifelse(toppar_sf$source == "fine", 2L, 1L)

# Välj vinnare per kluster
toppar_sf <- toppar_sf |>
  dplyr::group_by(grp) |>
  dplyr::arrange(dplyr::desc(prio), dplyr::desc(h_chm)) |>
  dplyr::slice(1) |>
  dplyr::ungroup() |>
  dplyr::select(geometry) |>
  dplyr::distinct(geometry, .keep_all = TRUE)

# Löpande ID
toppar_sf$treeID <- seq_len(nrow(toppar_sf))

# ====== 3d) Sätt topparnas CRS exakt till EPSG:3006 (och rensa Z/M) ======
toppar_sf <- sf::st_zm(toppar_sf, drop = TRUE, what = "ZM")
toppar_sf <- sf::st_transform(toppar_sf, target_epsg)  # säker reprojektion

# ====== 4) MCWS – kör via raster/sp med identisk CRS-objekt ======
# Skapa RasterLayer och SpatialPoints + tvinga samma CRS (WKT) på båda
chm_r <- raster::raster(chm_s)
crs_wkt <- sf::st_crs(target_epsg)$wkt           # WKT för EPSG:3006
raster::crs(chm_r) <- sp::CRS(SRS_string = crs_wkt)

treetops_sp <- methods::as(toppar_sf, "Spatial")
sp::proj4string(treetops_sp) <- sp::CRS(SRS_string = crs_wkt)

# Kör mcws
crowns_sp <- ForestTools::mcws(
  treetops  = treetops_sp,
  CHM       = chm_r,
  minHeight = 3,
  format    = "polygons"
)

# Tillbaka till sf
crowns <- sf::st_as_sf(crowns_sp)

# ====== 5) Efterbehandling ======
crowns <- crowns[ sf::st_area(crowns) > set_units(5, m^2), ]

# ====== 6) Spara ======
sf::st_write(
  crowns,
  file.path("C:/Path/For/Output.shp"),
  delete_layer = TRUE
)
sf::st_write(
  toppar_sf,
  file.path("C:/Path/For/Output.shp"),
  delete_layer = TRUE
)


# (valfritt) rensa temporärer
terra::tmpFiles(current = TRUE, remove = TRUE)

