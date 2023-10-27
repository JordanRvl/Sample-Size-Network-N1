library(data.table)

# Gather file names
folder = "data/data_sim/"
files = list.files(folder, pattern = "*.csv")

# Merging datasets from the different files
df_sim = list()
for (f in 1:length(files)) {
    df_sim[[f]] = read.csv(paste0(folder, files[f]), header=TRUE)
}
data = rbindlist(df_sim, fill=TRUE)

# Export
write.csv(data, "data/df_results.csv", row.names=FALSE)
print("MERGED")
