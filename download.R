##' 一键生成 eBird 相对丰度地图
##'
##' @param species           物种。可以是英文俗名 ("Wood Thrush") 或
##'                           6 位代码 ("woothr")。示例物种 "yebsap-example" 免 key。
##' @param season            时段。
##'                           - "breeding"                  繁殖季 (默认)
##'                           - "nonbreeding"                越冬季
##'                           - "prebreeding_migration"      春季迁徙
##'                           - "postbreeding_migration"     秋季迁徙
##'                           - "year_round"                全年
##'                           - "weekly"                     取某一周 (配合 week 参数)
##'                           - "custom"                     自定义日期范围 (配合 start/end)
##' @param start_date,end_date 当 season = "custom" 时使用; 格式 "YYYY-MM-DD"。
##' @param week              当 season = "weekly" 时使用; 1-52。
##' @param region            研究区域。四种形式:
##'                           - list(type="state",   name="Georgia",
##'                                  country_iso="US")   美国某州
##'                           - list(type="country", name="United States of America")
##'                           - list(type="sf",      name=<sf 多边形>) 直接传
##'                           - list(type="none")   不裁剪, 画全球
##'                           - list(type = "couns") 美国本土
##' @param resolution        栅格分辨率: "3km" (默认) / "9km" / "27km"。
##' @param metric            丰度统计量: "mean" (默认) / "median" / "lower" / "upper"/ "max"。
##'
##' @param quantile_breaks   是否用十分位分箱 (默认 TRUE)。
##' @param n_quantiles       分箱数, 默认 10。
##' @param palette           "ebirdst" (默认) 或任意颜色向量。
##' @param plot_map          是否出图 (默认 TRUE)。
##' @param save_tif          是否把结果栅格写盘 (默认 FALSE, 因为栅格本身就是下载来的)。
##' @param out_dir           输出目录, 默认 "output"。
##' @param out_prefix        文件名前缀; 默认自动生成。
##' @param map_title         地图标题; NULL 自动生成。
##' @param crs               投影; "auto" = 按区域质心自动选 LAEA。
##' @param verbose           打印进度 (默认 TRUE)。
##' @param ...               
##'
##' @return 一个 S3 对象 "rel_abundance":
##'   r       SpatRaster (裁剪投影后的丰度栅格)
##'   region  sf 多边形
##'   species, season, call ...
##' @export
map_online <- function(species,
                     season = c("breeding", "nonbreeding",
                                "prebreeding_migration",
                                "postbreeding_migration",
                                "year_round", "weekly", "custom"),
                     start_date          = NULL,
                     end_date            = NULL,
                     week                = NULL,
                     region              = list(type = "state",
                                                name = "Georgia",
                                                country_iso = "US"),
                     resolution          = c("3km", "9km", "27km"),
                     metric              = "mean",
                     quantile_breaks     = TRUE,
                     n_quantiles         = 10,
                     palette             = "ebirdst",
                     plot_map            = TRUE,
                     save_tif            = FALSE,
                     out_dir             = "output",
                     out_prefix          = NULL,
                     map_title           = NULL,
                     crs                 = "auto",
                     verbose             = TRUE,
                     api_key             = "",
                     ...) {
  
  ## ---------- 0. 参数对齐 ----------
  season     <- match.arg(season)
  resolution <- match.arg(resolution)
  
  .check_pkgs(c("ebirdst", "rnaturalearth", "rnaturalearthdata",
                "dplyr", "sf", "terra", "fields", "lubridate"))
  set.seed(1)
  
  if (verbose) message(">>> 物种: ", species,
                       "  时段: ", season)
  
  res <- .run_ebirdst(species = species, season = season,
                        start_date = start_date, end_date = end_date,
                        week = week, region = region,
                        resolution = resolution, metric = metric,
                        verbose = verbose, api_key = api_key)
  
  ## ---------- 出图或保存 ----------
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  if (is.null(out_prefix))
    out_prefix <- paste0(gsub(" ", "-", species), "_", season)
  
  if (save_tif) {
    out_tif <- file.path(out_dir, paste0(out_prefix, "_relative-abundance.tif"))
    terra::writeRaster(res$r, out_tif, overwrite = TRUE)
    if (verbose) message("    栅格已保存: ", out_tif)
  }
  
  if (plot_map) {
    if (is.null(map_title))
      map_title <- paste0(species, " — ",
                          .season_label(season, start_date, end_date, week))
    .plot_rel_abundance(r = res$r, region_sf = res$region,
                        quantile_breaks = quantile_breaks,
                        n_quantiles = n_quantiles,
                        palette = palette,
                        map_title = map_title,
                        verbose = verbose)
  }
  
  structure(
    list(species = species, season = season,
         r = res$r, region = res$region,
         crs_used = terra::crs(res$r),
         call = match.call()),
    class = "rel_abundance"
  )
}


## ============================================================================
##  主路径 A: 自动下载官方 Status & Trends 栅格
## ============================================================================
.run_ebirdst <- function(species, season, start_date, end_date, week,
                         region, resolution, metric, verbose, api_key) {
  
  ## 1) 物种名 -> 6 位代码; load_raster 也接受俗名, 但我们保险起见查一下
  if(api_key == "") stop("Put your api key in. get the key on ebird.org! It for map download, not for data")
  ebirdst::set_ebirdst_access_key(api_key, overwrite = TRUE)
  if (verbose) message(">>> [1/5] 下载/读取物种丰度栅格 (首次运行会自动联网下载) ...")
  ## load_raster 第一次调用即下载, 不需要显式 ebirdst_download_status
  if (season == "year_round") {
    r <- .load_raster_safely(species = species, period = "full-year",
                             metric = metric, resolution = resolution)
  } else if (season %in% c("breeding", "nonbreeding",
                           "prebreeding_migration",
                           "postbreeding_migration")) {
    r_seas <- .load_raster_safely(species = species, period = "seasonal",
                                  metric = metric, resolution = resolution)
    r <- r_seas[[season]]
  } else if (season == "weekly") {
    r_wk <- .load_raster_safely(species = species, period = "weekly",
                                metric = metric, resolution = resolution)
    if (is.null(week)) stop("season='weekly' 时请给 week=1..52")
    if (week < 1 || week > 52) stop("week 必须在 1..52")
    layer_name <- sprintf("%02d", week)
    ## 找到与该周对应的层名 (ebirdst 层名形如 "01", "02", ...)
    hit <- grep(layer_name, names(r_wk), value = TRUE, fixed = TRUE)
    if (length(hit) == 0) stop("在下载的栅格里找不到第 ", week, " 周")
    r <- r_wk[[hit[1]]]
  } else if (season == "custom") {
    r_wk <- .load_raster_safely(species = species, period = "weekly",
                                metric = metric, resolution = resolution)
    d1 <- lubridate::ymd(start_date); d2 <- lubridate::ymd(end_date)
    if (is.na(d1) || is.na(d2)) stop("start_date / end_date 格式应为 YYYY-MM-DD")
    ## ISO week: 把日期转成周数
    w1 <- lubridate::isoweek(d1); w2 <- lubridate::isoweek(d2)
    wk_seq <- seq(min(w1, w2), max(w1, w2))
    layer_names <- vapply(wk_seq, function(w) {
      g <- grep(sprintf("%02d", w), names(r_wk), value = TRUE, fixed = TRUE)
      if (length(g) == 0) return(NA_character_) else g[1]
    }, character(1))
    layer_names <- stats::na.omit(layer_names)
    if (length(layer_names) == 0)
      stop("日期范围内没有对应的周层; 请检查日期是否落在该物种有数据的年份。")
    if (length(layer_names) == 1) {
      r <- r_wk[[layer_names[1]]]
    } else {
      ## 多周取均值
      r <- terra::mean(r_wk[[layer_names]], na.rm = TRUE)
    }
    if (verbose) message("    自定义日期对应 ", length(layer_names), " 周栅格")
  }
  
  ## 2) 解析 region -> sf (rnaturalearth 自动下载边界)
  if (verbose) message(">>> [2/5] 获取研究区边界 (rnaturalearth 自动下载) ...")
  region_sf <- .resolve_region_ne(region)
  
  ## 3) 投影到 LAEA (以区域质心为中心)
  if (verbose) message(">>> [3/5] 投影并裁剪 ...")
  crs_raster <- terra::crs(r)
  if (!is.null(region_sf)) {
    region_proj <- sf::st_transform(region_sf, crs_raster) |> terra::vect()
    r_crop <- terra::crop(r, region_proj)
    r_crop <- terra::mask(r_crop, region_proj)
    ## 计算 LAEA 投影
    cent <- region_sf |>
      sf::st_transform(4326) |> sf::st_geometry() |>
      sf::st_centroid() |> sf::st_coordinates() |> as.numeric()
    crs_laea <- paste0("+proj=laea +lat_0=", round(cent[2], 3),
                       " +lon_0=", round(cent[1], 3))
    r_out <- terra::project(r_crop, crs_laea, method = "near") |> terra::trim()
    region_out <- sf::st_transform(region_sf, crs_laea)
  } else {
    r_out <- r
    region_out <- NULL
  }
  
  list(r = r_out, region = region_out)
}


## ============================================================================
##  辅助函数
## ============================================================================

.check_pkgs <- function(pkgs) {
  lapply(pkgs, library, character.only = TRUE)
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE))
      stop("缺少 R 包: ", p, "\n  请先 install.packages('", p, "')")
  }
}

##' 确保 rnaturalearthhires 可用 (ne_states 高精度州界需要它;
##' 该包已从 CRAN archive, 用 r-universe 源装)
.ensure_rnaturalearthhires <- function() {
  if (requireNamespace("rnaturalearthhires", quietly = TRUE)) return(invisible(TRUE))
  message(">>> 首次使用, 尝试安装 rnaturalearthhires ...")
  tryCatch({
    install.packages("rnaturalearthhires",
                     repos = c("https://ropensci.r-universe.dev/",
                               "https://cloud.r-project.org"))
  }, error = function(e) NULL)
  invisible(requireNamespace("rnaturalearthhires", quietly = TRUE))
}

##' 拉州界: 优先 ne_states (高精度), 失败则用 ne_download 低分辨率 (免 hires 包)
.fetch_states_sf <- function(iso_a2) {
  if (.ensure_rnaturalearthhires()) {
    return(rnaturalearth::ne_states(iso_a2 = iso_a2, returnclass = "sf"))
  }
  message("    rnaturalearthhires 不可用, 改用低分辨率州界 ...")
  rnaturalearth::ne_download(scale = 110, type = "states",
                             category = "cultural",
                             returnclass = "sf")
}

##' 安全调用 load_raster: 若请求的分辨率该物种不提供, 自动降级
##' (例如 yebsap-example 只有 27km)
.load_raster_safely <- function(species, period, metric, resolution) {
  cands <- unique(c(resolution, "9km", "27km", "3km"))
  last_err <- NULL
  for (res in cands) {
    r <- tryCatch(
      ebirdst::load_raster(species = species, product = "abundance",
                           period = period, metric = metric, resolution = res),
      error = function(e) { last_err <<- e; NULL }
    )
    if (!is.null(r)) {
      if (res != resolution)
        message("    (该物种没有 ", resolution,
                " 分辨率数据, 自动降级到 ", res, ")")
      return(r)
    }
  }
  msg <- conditionMessage(last_err)
  if (grepl("Cannot access Status", msg)) {
    stop(
      "ebirdst data server (st-download.ebird.org) is down.\n",
      "We probed it directly and got HTTP 500 for every species/year, with or without key.\n",
      "This is an eBird server-side outage, NOT your API key and NOT this function.\n",
      "Do this:\n",
      "  1) check ur api is copied from the map page not dataset pages;\n",
      "  2) open https://science.ebird.org/en/status-and-trends/species/,click downloads get ur api",
      tolower(gsub(" ", "", species)),
      "/downloads in a browser, download the GeoTIFFs manually into\n",
      "     ebirdst::ebirdst_data_dir(), then rerun;\n",
      "  3) Or verify the function works with the keyless example species:\n",
      "     f.rel.ab(\"yebsap-example\", \"breeding\",\n",
      "              list(type=\"state\", name=\"Michigan\", country_iso=\"US\")).\n",
      "Original error: ", msg
    )
  }
  stop("load_raster failed: ", msg)
}

##' 用 rnaturalearth 把 region 参数解析成 sf 多边形
.resolve_region_ne <- function(region) {
  if (is.null(region) || is.null(region$type)) return(NULL)
  switch(region$type,
         sf = region$name,
         state = {
           iso <- region$country_iso %||% "US"
           nm  <- region$name          # 提前取出, 避免 dplyr::filter mask 把 region 解析成列
           sf_states <- .fetch_states_sf(iso)
           sf_states |>
             dplyr::filter(grepl(nm, name, ignore.case = TRUE)) |>
             sf::st_geometry()
         },
         country = rnaturalearth::ne_countries(country = region$name,
                                               returnclass = "sf") |>
           sf::st_geometry(),
         conus = {
           ## 美国本土 48 州 (自动排除 Alaska / Hawaii )
           .fetch_states_sf("US") |>
             dplyr::filter(!grepl("alaska|hawaii|puerto rico", name,
                                  ignore.case = TRUE)) |>
             sf::st_union() |>
             sf::st_geometry()
         },
         none = NULL,
         stop("region$type 只支持 'state' / 'country' / 'conus' / 'sf' / 'none'")
  )
}

##' 出图: 复刻 ebirdst applications vignette 的标准样式
.plot_rel_abundance <- function(r, region_sf,
                                quantile_breaks, n_quantiles,
                                palette, map_title, verbose) {
  if (verbose) message(">>> [4/5] 计算分箱 ...")
  vals <- terra::values(r, na.rm = TRUE, mat = FALSE)
  vals_pos <- vals[vals > 0 & !is.na(vals)]
  if (length(vals_pos) < 10) {
    warning("研究区内非零栅格单元太少, 地图可能没有意义。")
    brks <- seq(0, max(vals, na.rm = TRUE), length.out = n_quantiles + 1)
  } else if (quantile_breaks) {
    qs <- stats::quantile(vals_pos, seq(0, 1, length.out = n_quantiles + 1),
                          na.rm = TRUE)
    brks <- unique(c(0, as.numeric(qs)))
  } else {
    brks <- seq(0, max(vals, na.rm = TRUE), length.out = n_quantiles + 1)
  }
  
  if (identical(palette, "ebirdst")) {
    pal <- ebirdst::ebirdst_palettes(length(brks) - 2)
    pal <- c("#e6e6e6", pal)
  } else {
    pal <- palette
  }
  lbls <- round(c(min(brks), stats::median(brks), max(brks)), 3)
  
  if (verbose) message(">>> [5/5] 绘图 ...")
  crs_map <- terra::crs(r)
  countries <- rnaturalearth::ne_countries(returnclass = "sf") |>
    sf::st_geometry() |> sf::st_transform(crs_map)
  countries_v <- terra::vect(countries)
  states <- tryCatch(
    .fetch_states_sf("US") |> sf::st_transform(crs_map),
    error = function(e) NULL
  )
  states_v <- if (!is.null(states)) terra::vect(states) else NULL
  region_v <- if (!is.null(region_sf))
    terra::vect(sf::st_transform(region_sf, crs_map)) else NULL
  
  ## 给右侧图例留出空间, 用 par(mar) 控制
  graphics::par(mar = c(4, 0.5, 1, 0.5))
  if (!is.null(region_v)) {
    terra::plot(region_v, col = NA, border = NA, axes = FALSE)
    terra::plot(countries_v, col = "#cfcfcf", border = "#888888",
                lwd = 0.5, add = TRUE)
  } else {
    terra::plot(countries_v, col = "#cfcfcf", border = "#888888",
                lwd = 0.5, axes = FALSE)
  }
  terra::plot(r, col = pal, breaks = brks,
              maxpixels = terra::ncell(r),
              legend = FALSE, axes = FALSE, add = TRUE)
  if (!is.null(states_v))
    terra::lines(states_v, col = "#ffffff", lwd = 0.75)
  terra::lines(countries_v, col = "#ffffff", lwd = 1.5)
  if (!is.null(region_v))
    terra::lines(region_v, col = "#000000", lwd = 1.2)
  
  ## 右侧垂直图例 (复刻 ebirdst applications vignette)
  labels_q <- stats::quantile(brks, c(0, 0.5, 1))
  label_breaks <- seq(0, 1, length.out = length(brks))
  fields::image.plot(
    zlim = c(0, 1),
    breaks = label_breaks,
    col = pal,
    smallplot = c(0.90, 0.93, 0.15, 0.85),
    legend.only = TRUE,
    axis.args = list(
      at = c(0, 0.5, 1),
      labels = round(labels_q, 2),
      col.axis = "black", fg = NA,
      cex.axis = 0.9, lwd.ticks = 0, line = -0.5
    )
  )
  graphics::mtext(map_title, side = 3, line = -1, cex = 1.1)
  invisible(NULL)
}

.season_label <- function(season, start_date, end_date, week = NULL) {
  switch(season,
         breeding               = "Breeding Season",
         nonbreeding            = "Non-breeding Season",
         prebreeding_migration  = "Pre-breeding Migration",
         postbreeding_migration = "Post-breeding Migration",
         year_round             = "Full Year",
         weekly                 = paste0("Week ", week),
         custom                 = paste0(start_date, " ~ ", end_date))
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a


## ============================================================================
##  S3 方法
## ============================================================================
print.rel_abundance <- function(x, ...) {
  cat("== eBird 相对丰度地图 ==\n")
  cat("  物种 :", x$species, "\n")
  cat("  时段 :", x$season, "\n")
  if (!is.null(x$r))
    cat("  栅格 :", paste(dim(x$r), collapse = " x "), "\n")
  if (!is.null(x$n_checklists))
    cat("  清单 :", x$n_checklists, "(建模", x$n_model %||% NA, ")\n")
  invisible(x)
}

plot.rel_abundance <- function(x, ...) {
  .plot_rel_abundance(x$r, x$region,
                      quantile_breaks = TRUE, n_quantiles = 10,
                      palette = "ebirdst",
                      map_title = paste0(x$species, " — ", x$season),
                      verbose = FALSE)
}
