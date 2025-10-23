############################################################
# 万博：Visitors × Temperature
# 第1〜4回の範囲：可視化 → ADF/差分 → ACF/PACF → VAR → 因果/IRF/FEVD → 診断/予測
############################################################

# 必要パッケージ（未導入なら install.packages(...) で導入）
# install.packages(c("lubridate","tseries","urca","vars","forecast","PerformanceAnalytics","ggplot2"))

library(lubridate)
library(tseries)
library(urca)
library(vars)
library(forecast)
library(PerformanceAnalytics)
library(ggplot2)

# ---- 1) データ読み込み & 整形（第1回：時系列の準備・可視化） ----
csv_path <- "/root/workspace/ToDM/dataset/expo_visitor_data.csv"

raw <- read.csv(csv_path, stringsAsFactors = FALSE)
names(raw) <- c("date", "month", "day", "weekday", "visitor_count", "weather", "temperature")
stopifnot(all(c("date","visitor_count","temperature") %in% names(raw)))

df <- data.frame(
  Date = as.Date(raw$date),
  Visitors = as.numeric(raw$visitor_count),
  Temperature = as.numeric(raw$temperature)
)
df <- df[order(df$Date), ]
df <- df[!is.na(df$Date) & !is.na(df$Visitors) & !is.na(df$Temperature), ]

# 基本可視化（第1回）
dir.create("/root/workspace/ToDM/fig", recursive = TRUE, showWarnings = FALSE)
png("/root/workspace/ToDM/fig/visitors.png", 800, 600)
plot(df$Date, df$Visitors, type="l", main="Daily Visitors", xlab="Date", ylab="Visitors")
dev.off()
png("/root/workspace/ToDM/fig/temperature.png", 800, 600)
plot(df$Date, df$Temperature, type="l", main="Daily Temperature", xlab="Date", ylab="Temperature")
dev.off()

# ---- 2) 変換・定常性（第2〜3回：ADF と差分） ----
df$Visitors_log <- log1p(df$Visitors)  # 分散安定化の一例（第2回の変換の考え方）
df$Temp <- df$Temperature

adf_p <- function(x) suppressWarnings(tseries::adf.test(x, k = trunc((length(x)-1)^(1/3))))$p.value

make_stationary <- function(vec, max_diff = 2, name="series") {
  lvl <- vec; d <- 0; p <- adf_p(lvl)
  while (p > 0.05 && d < max_diff) {
    lvl <- diff(lvl); d <- d + 1; p <- adf_p(lvl)
  }
  message(sprintf("[%s] differenced %d time(s); ADF p=%.4f", name, d, p))
  list(x = lvl, d = d, pval = p)
}

vis_res <- make_stationary(df$Visitors_log, name="Visitors_log")  # 第2・3回
tmp_res <- make_stationary(df$Temp,         name="Temperature")

# 長さ合わせ
align_series <- function(a, b) { n <- min(length(a), length(b)); cbind(a = tail(a, n), b = tail(b, n)) }
Z <- align_series(vis_res$x, tmp_res$x)
colnames(Z) <- c("Visitors_s", "Temp_s")
end_dates <- tail(df$Date, nrow(as.data.frame(Z)))
ts_df <- data.frame(Date = end_dates, Visitors_s = as.numeric(Z[,1]), Temp_s = as.numeric(Z[,2]))

# ACF/PACF（第3回：構造把握とラグの目安）
png("/root/workspace/ToDM/fig/acf_pacf_visitors_s.png", 1000, 800)
par(mfrow=c(2,1)); acf(ts_df$Visitors_s, main="ACF Visitors_s"); pacf(ts_df$Visitors_s, main="PACF Visitors_s"); par(mfrow=c(1,1))
dev.off()
png("/root/workspace/ToDM/fig/acf_pacf_temp_s.png", 1000, 800)
par(mfrow=c(2,1)); acf(ts_df$Temp_s, main="ACF Temp_s"); pacf(ts_df$Temp_s, main="PACF Temp_s"); par(mfrow=c(1,1))
dev.off()

# ---- 3) VAR 次数選択・推定（第3回：情報量規準／第4回：VAR）----
sel <- VARselect(ts_df[, c("Visitors_s","Temp_s")], lag.max = 14, type = "const")
print(sel$criteria)
p_opt <- if(!is.na(sel$selection["AIC(n)"])) sel$selection["AIC(n)"] else sel$selection[1]
var_fit <- VAR(ts_df[, c("Visitors_s","Temp_s")], p = p_opt, type = "const")
print(summary(var_fit))

# ---- 4) 診断（第3回：診断）＋（第4回：多変量残差の検定“まで”）----
cat("---- Serial correlation (Portmanteau) ----\n")
print(serial.test(var_fit, lags.pt = 12, type = "PT.asymptotic"))

cat("---- ARCH effects (multivariate ARCH LM) ----\n")
print(arch.test(var_fit, lags.multi = 12))  # ここは“検定”まで（推定はしない）

cat("---- Normality test (Jarque-Bera) ----\n")
print(normality.test(var_fit))

# ---- 5) グレンジャー因果・IRF・FEVD（第4回）----
cat("---- Granger Causality ----\n")
print(causality(var_fit, cause = "Temp_s"))
print(causality(var_fit, cause = "Visitors_s"))

cat("---- IRF ----\n")
irf_temp_to_vis <- irf(var_fit, impulse = "Temp_s", response = c("Visitors_s","Temp_s"), n.ahead = 8, boot = TRUE)
png("/root/workspace/ToDM/fig/irf_temp_to_visitors.png", 1000, 800); plot(irf_temp_to_vis); dev.off()
irf_vis_to_temp <- irf(var_fit, impulse = "Visitors_s", response = c("Visitors_s","Temp_s"), n.ahead = 8, boot = TRUE)
png("/root/workspace/ToDM/fig/irf_visitors_to_temp.png", 1000, 800); plot(irf_vis_to_temp); dev.off()

cat("---- FEVD ----\n")
fevd_res <- fevd(var_fit, n.ahead = 8)
png("/root/workspace/ToDM/fig/fevd.png", 1000, 800); plot(fevd_res); dev.off()

# ---- 6) 短期予測（第3回：予測）----
h <- 7
var_fc <- predict(var_fit, n.ahead = h, ci = 0.95)
print(var_fc)

vis_fc <- data.frame(
  step  = 1:h,
  fc    = var_fc$fcst$Visitors_s[, "fcst"],
  lower = var_fc$fcst$Visitors_s[, "lower"],
  upper = var_fc$fcst$Visitors_s[, "upper"]
)
p_fc <- ggplot(vis_fc, aes(step, fc)) +
  geom_line() + geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.25) +
  labs(title = sprintf("VAR(%d) %d-step ahead forecast (Visitors_s)", p_opt, h),
       x = "Steps ahead (days)", y = "Forecast") + theme_minimal()
print(p_fc)
ggsave("/root/workspace/ToDM/fig/var_forecast_7days.png", p_fc, width = 10, height = 6, dpi = 300)

# ---- 7) 簡易サマリ ----
resid_mat <- residuals(var_fit)
cat("---- Model Summary ----\n")
cat("VAR order (AIC):", p_opt, "\n")
cat("Obs:", nrow(ts_df), "\n")
cat("Residual sd (Visitors_s):", sd(resid_mat[, "Visitors_s"]), "\n")
cat("Residual sd (Temp_s):",     sd(resid_mat[, "Temp_s"]), "\n")