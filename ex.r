## ============================================================================
##  ex.r —— 两个出图函数的用法示例
## ----------------------------------------------------------------------------
##  本工程有两个并列的出图函数, 各管一条路线:
##
##  【函数 A: map_local】  —— 在 custom.r 里
##    用你自己下载的 eBird 基本数据集 (EBD/SED 两个 .txt), 本地从零训练
##    hurdle 模型再出相对丰度图。完全自主、不依赖 eBird 发布的成品产品;
##    环境栅格 (高程/土地覆盖) 首次运行自动联网下载并缓存。
##
##  【函数 B: map_ebirdst】 —— 你单独维护的那个函数 (负责在线下载官方成品栅格)
##    直接下载 eBird Status & Trends 已发布的相对丰度栅格, 快, 需要 eBird key。
##    注意: 它的函数名在你另一个文件里, 下面示例里用 map_ebirdst 占位,
##          请把 "map_ebirdst" 改成你实际保存的那个函数名 (参数照抄即可)。
##
##  建议在 RStudio 里打开本工程的 .Rproj, 工作目录即工程根, 相对路径才生效。
## ============================================================================

## ---- 加载函数 A (本地建模) --------------------------------------------------
source("custom.r")


## ============================================================================
##  函数 A: map_local() —— 用本地 EBD/SED 从零建模
##  准备: 把 EBD/SED 两个纯文本放进 data/ (子文件夹也行, 自动递归找)。
## ============================================================================

## ---- A1【最稳, 第一次先跑这个】---------------------------------------------
## region=auto 按数据范围定研究区; season=year_round 不挑月份, 最不容易失败。
## 不确定数据里有哪些鸟: 故意填错名字, 报错会列出数据里最常见的鸟, 照抄即可。
a1 <- map_local(
  species = "Northern Cardinal",        # 通用名; 也可填学名 "Cardinalis cardinalis"
  season  = "year_round",               # 全年, 最稳
  region  = list(type = "auto")         # 研究区 = 数据覆盖范围
)

## ---- A2: 指定季节 (要和你数据的月份对得上, 对不上会自动放宽到全年并提示) ----
a2 <- map_local(
  species = "Northern Cardinal",
  season  = "prebreeding_migration",    # breeding/nonbreeding/pre/post/year_round/weekly/custom
  region  = list(type = "auto")
)

## ---- A3: 一年中的第 N 周 ---------------------------------------------------
a3 <- map_local(
  species = "Northern Cardinal",
  season  = "weekly",
  week    = 12,                         # 第 12 周 (约 3 月中下旬)
  region  = list(type = "auto")
)

## ---- A4: 完全自定义日期区间 -------------------------------------------------
a4 <- map_local(
  species     = "Northern Cardinal",
  season      = "custom",
  start_date  = "2025-03-05",
  end_date    = "2025-03-20",
  region      = list(type = "auto")
)

## ---- A5: 指定行政区域 (数据须落在区内, 否则自动放宽到数据范围) --------------
a5a <- map_local(species = "Northern Cardinal", season = "year_round",
                 region = list(type = "state", name = "Alabama", country_iso = "US"))
a5b <- map_local(species = "Northern Cardinal", season = "year_round",
                 region = list(type = "country", name = "United States of America"))
a5c <- map_local(species = "Northern Cardinal", season = "year_round",
                 region = list(type = "conus"))        # 美国本土 48 州
## 也可直接传一个经纬度边界框 (sf 对象)
my_bbox <- sf::st_bbox(c(xmin = -87.8, ymin = 33.0, xmax = -87.0, ymax = 33.5),
                       crs = sf::st_crs(4326))
a5d <- map_local(species = "Northern Cardinal", season = "year_round",
                 region = sf::st_as_sfc(my_bbox))

## ---- A6: 学名输入 -----------------------------------------------------------
a6 <- map_local(species = "Cardinalis cardinalis", season = "year_round",
                region = list(type = "auto"))

## ---- A7: 控制精细度 / 建模细节 ---------------------------------------------
a7 <- map_local(
  species        = "Northern Cardinal",
  season         = "prebreeding_migration",
  region         = list(type = "auto"),
  grid_res       = 3000,    # 预测网格边长(米); 小=更细更噪, 大=更平滑
  elev_z         = 7,       # 高程精度等级
  lc_agg         = 90,      # 土地覆盖聚合(米), 大=快而平滑
  sample_radius  = 1500,    # 点/网格周围环境统计半径(米)
  n_trees        = 500,     # 随机森林树数
  smooth         = TRUE      # 小样本空间去噪(默认开); FALSE 完全复刻教程
)

## ---- A8: 数据不在默认 data/ 时, 手动指定路径 -------------------------------
a8 <- map_local(
  species  = "Northern Cardinal",
  season   = "year_round",
  region   = list(type = "auto"),
  data_dir = "D:/my_ebird_data",
  out_dir  = "my_maps"
)


## ============================================================================
##  函数 B: map_ebirdst() —— 在线下载官方成品栅格 (快, 需 eBird key)
##  !!! 下面的 map_ebirdst 是占位函数名: 请改成你单独维护的那个函数的名字 !!!
##  它的参数与原来 f.rel.ab(method="ebirdst", ...) 一致。
## ============================================================================

## 免 key 示例物种, 用来验证这条链路是否正常
b_demo <- map_ebirdst(
  species    = "yebsap-example",
  season     = "breeding",
  region     = list(type = "state", name = "Michigan", country_iso = "US"),
  metric     = "median",        # median / mean / upper / lower
  resolution = "3km"            # 3km 或 27km
)

## 真实物种示例
b1 <- map_ebirdst(
  species    = "Wood Thrush",
  season     = "breeding",
  region     = list(type = "state", name = "Georgia", country_iso = "US"),
  metric     = "median",
  resolution = "3km"
)

## 周度 / 自定义日期
b2 <- map_ebirdst(species = "Wood Thrush", season = "weekly", week = 24,
                  region = list(type = "country", name = "United States of America"),
                  metric = "upper")


## ============================================================================
##  返回结果与输出文件 (两个函数返回同一种 S3 对象)
## ============================================================================
## 地图默认立即弹出, 并在 out_dir(默认 "output/") 保存:
##   * PDF 地图
##   * <前缀>_abundance.tif —— 最终相对丰度
##   * <前缀>_layers.tif    —— 四层: in_range / encounter_rate / count / abundance
print(a1)
plot(a1)
terra::plot(a1$r_layers[["encounter_rate"]])   # 只看遭遇率层
a1$diagnostics       # 模型诊断(阈值/测试误差等)
a1$season_used       # 实际使用的季节(被自动放宽时会和你填的不同)
