sink("Result_taxEffort.txt", split = TRUE)

library(QCA); library(SetMethods); library(dplyr); library(readxl)

data <- read_excel("grb_panel.xlsx")

q   <- function(x, p) as.numeric(quantile(x, p, na.rm = TRUE))
cal <- function(x) calibrate(x, type = "fuzzy",
                             thresholds = c(q(x, .05), q(x, .5), q(x, .95)))

cd <- data %>% transmute(
  SIDO    = sido,
  YEAR    = year,
  OUT     = cal(taxEffort),             # tax collection rate
  FEM_CS  = cal(per_emp_women),         # share of female civil servants
  LIB_MAY = ifelse(liberal_may == 1, 1, 0), # liberal mayor
  FIS_AUT = cal(fisAuto),               # fiscal autonomy
  POP     = cal(lnpoptotal),            # log population
  CAP_REG = ifelse(cap_reg == 1, 1, 0)  # capital region
)

L  <- c("FEM_CS", "LIB_MAY", "FIS_AUT", "POP", "CAP_REG")
cd <- as.data.frame(cd[complete.cases(cd[c("OUT", L)]), ])
rownames(cd) <- paste(cd$SIDO, cd$YEAR)
cat("N =", nrow(cd), " units =", length(unique(cd$SIDO)),
    " years =", length(unique(cd$YEAR)), "\n")

con  <- function(x, y) (sum(pmin(x, y)) + 1e-10) / (sum(x) + 1e-10)
cvg  <- function(x, y) (sum(pmin(x, y)) + 1e-10) / (sum(y) + 1e-10)
byg  <- function(f, x, y, g) { k <- sort(unique(g))
setNames(sapply(k, function(i) f(x[g == i], y[g == i])), k) }
adjd <- function(z, n) {
  z <- (z + 1e-10) / (sum(z) + 1e-10)
  sqrt(sum((z - 1/length(z))^2)) / sqrt(n / (n^2 + 3*n + 2))
}
nyr <- length(unique(cd$YEAR)); nun <- length(unique(cd$SIDO))


pres <- as.data.frame(lapply(cd[L], as.numeric))
absn <- as.data.frame(lapply(cd[L], function(x) 1 - as.numeric(x)))
names(absn) <- paste0("~", L)

for (y in list(OUT = cd$OUT, notOUT = 1 - cd$OUT)) {
  print(pof(pres, y, relation = "necessity"))
  print(pof(absn, y, relation = "necessity"))
}

for (nm in c("OUT", "~OUT")) {
  Y <- if (nm == "OUT") cd$OUT else 1 - cd$OUT
  tab <- t(sapply(c(L, paste0("~", L)), function(v) {
    x <- cd[[sub("^~", "", v)]]; if (startsWith(v, "~")) x <- 1 - x
    be <- byg(cvg, x, Y, cd$YEAR)
    c(POCONS = cvg(x, Y), be, BE_adjd = adjd(be, nyr))
  }))
  cat("\n== Necessity for", nm, ": POCONS, BECONS by year, distance (cut .90) ==\n")
  print(round(tab, 3))
}


# Sufficiency
tt <- truthTable(cd, outcome = "OUT", conditions = L,
                 incl.cut = .8, pri.cut = .6, n.cut = 3,
                 show.cases = FALSE, sort.by = "incl")
print(tt)

sol <- minimize(tt, include = "?", dir.exp = setNames(rep(1, length(L)), L),
                details = TRUE)

cat("\n=== COMPLEX solution ===\n");      print(minimize(tt, details = TRUE))
cat("\n=== PARSIMONIOUS solution ===\n"); print(minimize(tt, include = "?", details = TRUE))
cat("\n=== INTERMEDIATE solution ===\n"); print(sol)


# Cluster diagnostics
mods   <- unlist(lapply(names(sol$i.sol), function(cp)
  paste0(tolower(cp), "i", seq_along(sol$i.sol[[cp]]$solution))))
paths  <- unlist(lapply(names(sol$i.sol), function(cp) sol$i.sol[[cp]]$solution),
                 recursive = FALSE)
mods   <- mods[!duplicated(lapply(paths, sort))]

for (s in mods) {
  pd <- pimdata(results = sol, outcome = "OUT", sol = s)
  Y  <- cd$OUT
  cat("\n-- solution", s, "\n--")
  
  cat("\n-- Pooled scores and adjusted distances (cut .1 / .2) --\n")
  print(round(t(sapply(pd[-ncol(pd)], function(x) c(
    POCONS  = con(x, Y),
    POCOV   = cvg(x, Y),
    BE_adjd = adjd(byg(con, x, Y, cd$YEAR), nyr),
    WI_adjd = adjd(byg(con, x, Y, cd$SIDO), nun)))), 3))
  
  cat("\n-- BECONS by year --\n")
  print(round(t(sapply(pd[-ncol(pd)], function(x) byg(con, x, Y, cd$YEAR))), 3))
  cat("\n-- BECOV by year --\n")
  print(round(t(sapply(pd[-ncol(pd)], function(x) byg(cvg, x, Y, cd$YEAR))), 3))
  
  X  <- pd$solution_formula
  wc <- byg(con, X, Y, cd$SIDO); wv <- byg(cvg, X, Y, cd$SIDO)
  cat("\n-- WICONS distribution, solution as a whole (n =", nun, ") --\n")
  print(round(quantile(wc, c(0, .05, .25, .5, .75, .95, 1)), 3))
  cat(sprintf("  =1.00: %d | .80-1.00: %d | <.80: %d | <.50: %d   WICOV median = %.3f\n",
              sum(wc >= .9999), sum(wc >= .8 & wc < .9999), sum(wc < .8),
              sum(wc < .5), median(wv)))
  cat("\n-- Least consistent units (theory violations; WICOV = empirical relevance) --\n")
  o <- order(wc)[1:10]
  print(data.frame(unit = sort(unique(cd$SIDO))[o],
                   WICONS = round(wc[o], 3), WICOV = round(wv[o], 3)), row.names = FALSE)
}


# Ten cases for each sufficient configuration
for (s in mods) {
  typ <- smmr(results = sol, outcome = "OUT", sol = s, match = FALSE, cases = 1)[[1]]$results
  for (p in unique(as.character(typ$Term))) {
    tp <- typ[as.character(typ$Term) == p, ]
    cat("\n== ", p, "  (typical n = ", nrow(tp), " case-years) ==\n", sep = "")
    print(head(tp[order(-tp$TermMemb), c("Case", "TermMemb", "Outcome")], 10),
          row.names = FALSE)
  }
}

sink()