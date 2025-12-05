# Requires pacman, data.table, tidyverse, readxl, magrittr
pacman::p_load(data.table, tidyverse, readxl, magrittr);

# This script is just an example of how to run chisq.test on the metadata with simulated p-values.
# Adjust the path and the number of iterations below as you will.

# Path to metadata file
path <- '~/Downloads/New lifestyle survey fields.xlsx';

# Set number of MonteCarlo simulations for p-value
num.iter <- 100000;

# Read as excel file as tibble and convert to data.table
dt <- as.data.table(readxl::read_xlsx(path));

# Convert all character columns to factors
cols <- dt[, lapply(.SD, is.character)] %>% unlist() %>% which() %>% names();
dt[, (cols) := lapply(.SD, factor), .SDcols = cols];

# Create an empty list of chisq tests of Ethnic_1 and bread_sum_edited (comparison named ethnic.bread)
chisq.res <- list(ethnic.bread = list(analytic = NULL, simulate = NULL));

# Fill list with analytic chisq.test result.
# This should generate a warning because > 20% of data has expected count < 5.
chisq.res$ethnic.bread$analytic <- chisq.test(dt$Ethnic_1, dt$bread_sum_edited);

# Fill list with simulated chisq.test result.
# Default number of iterations (B=2000) is not enough.
# Minimum simulated p-value >= 1/(num_iter + 1), so the greater num_iter, the greater detection of very small p-values.
chisq.res$ethnic.bread$simulate <- chisq.test(dt$Ethnic_1, dt$bread_sum_edited, simulate.p.value = TRUE, B = num.iter);

# Show results
chisq.res$ethnic.bread$analytic
chisq.res$ethnic.bread$simulate
