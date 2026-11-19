# Example gallery for the plot/summary methods of "mlmridge" objects (multivariate ridge).
# Each method delegates to the univariate "ridge" method for one response via as.ridge().
# Writes PNGs plus an HTML page (code next to each plot) to `outdir`.
#
# Run from the package root:  Rscript dev/mlmridge-plot-examples.R [outdir]

suppressMessages(devtools::load_all(".", quiet = TRUE))
args <- commandArgs(TRUE)
outdir <- if (length(args)) args[1] else file.path(tempdir(), "mlmridge-examples")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

data(schooldata, package = "heplots")
lambda <- c(0, 0.1, 0.5, 1, 2, 5)
sridge <- ridge(cbind(reading, mathematics, selfesteem) ~
                education + occupation + visit + counseling + teacher,
                data = schooldata, lambda = lambda)

examples <- list(
  list(id = "traceplot", w = 1200, h = 450,
       title = "traceplot(): by default, all responses, one panel each",
       code = 'traceplot(sridge)'),
  list(id = "traceplot-df", w = 1200, h = 450,
       title = "traceplot(X = \"df\"): same, against effective degrees of freedom",
       code = 'traceplot(sridge, X = "df")'),
  list(id = "plot", w = 1200, h = 450,
       title = "plot(): bivariate ridge trace for occupation and visit, each response",
       code = 'op <- par(mfrow = c(1, 3))
for (r in c("reading", "mathematics", "selfesteem"))
  plot(sridge, response = r, variables = c("occupation", "visit"),
       radius = 0.5, labels = lambda)
par(op)'),
  list(id = "cross-one", w = 1200, h = 450,
       title = "plot(predictor = ): coefficient ellipses across two responses (reading, mathematics)",
       code = 'op <- par(mfrow = c(1, 3))
for (p in c("occupation", "visit", "counseling"))
  plot(sridge, predictor = p, responses = c("reading", "mathematics"), fill = TRUE)
par(op)'),
  list(id = "cross-all", w = 700, h = 650,
       title = "plot(predictor = 1:5): all predictors, at lambda = 0 (black) and 5 (red)",
       code = 'plot(sridge, predictor = 1:5, which.lambda = c(1, 6), col = c("black", "red"))'),
  list(id = "pairs", w = 800, h = 800,
       title = "pairs(): all pairs of predictors, response = \"reading\"",
       code = 'pairs(sridge, response = "reading", radius = 0.5)'),
  list(id = "biplot", w = 1200, h = 450,
       title = "pca() + plot(): the two smallest PCA dimensions, each response",
       code = 'op <- par(mfrow = c(1, 3))
for (r in c("reading", "mathematics", "selfesteem")) {
  pr <- pca(sridge, response = r)
  plot(pr, radius = 0.5, labels = lambda, main = r)
}
par(op)'),
  list(id = "vif", w = 700, h = 500,
       title = "vif(): the same for every response (depends only on X)",
       code = 'v <- vif(sridge)
plot(v)'),
  list(id = "precision", w = 1200, h = 450,
       title = "precision(): variance vs. shrinkage, each response",
       code = 'op <- par(mfrow = c(1, 3))
for (r in c("reading", "mathematics", "selfesteem")) {
  pdat <- precision(sridge, response = r)
  plot(pdat, main = r)
}
par(op)')
)

html <- c('<!doctype html><html><head><meta charset="utf-8"><title>mlmridge plot examples</title>',
  '<style>body{font-family:sans-serif;max-width:1250px;margin:auto;padding:1em}',
  'pre{background:#f4f4f4;padding:.7em;border-radius:5px} img{max-width:100%;border:1px solid #ddd}',
  'h2{font-size:1.1em;margin-top:2.2em}</style></head><body>',
  '<h1>Plot methods for <code>"mlmridge"</code> objects</h1>',
  '<p>Data: <code>heplots::schooldata</code>. Setup:</p>',
  '<pre>sridge &lt;- ridge(cbind(reading, mathematics, selfesteem) ~\n                education + occupation + visit + counseling + teacher,\n                data = schooldata, lambda = c(0, 0.1, 0.5, 1, 2, 5))</pre>')

for (ex in examples) {
  f <- paste0(ex$id, ".png")
  png(file.path(outdir, f), width = ex$w, height = ex$h, res = 110)
  msgs <- character()
  withCallingHandlers(eval(parse(text = ex$code)),
                      message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
  dev.off()
  html <- c(html, sprintf("<h2>%s</h2>", ex$title),
            sprintf("<pre>%s</pre>", gsub("<", "&lt;", ex$code)),
            if (length(msgs)) sprintf("<p><i>message:</i> %s</p>", paste(trimws(msgs), collapse = "; ")),
            sprintf('<img src="%s">', f))
}

# the default-response message, and print/vif output, as text
txt <- capture.output({
  cat("> plot(sridge)   # no response = : uses the first response, with a message\n")
  withCallingHandlers({ pdf(NULL); plot(sridge); dev.off() },
                      message = function(m) { cat(conditionMessage(m)); invokeRestart("muffleMessage") })
  cat("\n> vif(sridge)\n"); print(vif(sridge))
})
html <- c(html, "<h2>Text output</h2>", sprintf("<pre>%s</pre>", paste(txt, collapse = "\n")),
          "</body></html>")
writeLines(html, file.path(outdir, "index.html"))
cat("wrote", length(examples), "plots and index.html to", outdir, "\n")
