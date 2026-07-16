packages <- c(
    "ggplot2",
    "dplyr",
    "patchwork",
    "ggExtra",
    "VennDiagram"
)
for (pkg in packages) {
    suppressPackageStartupMessages(
        suppressWarnings(
            library(
                pkg,
                character.only = TRUE,
                quietly = TRUE,
                warn.conflicts = FALSE
            )
        )
    )
}
source("../../utils/r_themes.r")

lm_results_file_path <- file.path(
    "../results/all_features_beta_df.parquet"
)
lm_coeff_df <- arrow::read_parquet(lm_results_file_path)
# shuffle the row order for plotting purposes
lm_coeff_df <- lm_coeff_df %>%
    dplyr::mutate(
        row_id = 1:nrow(lm_coeff_df)
    ) %>%
    dplyr::arrange(dplyr::desc(row_id)) %>%
    dplyr::select(-row_id)

lm_coeff_df$log10p_value <- -log10(lm_coeff_df$p_value)
# remove the const from the variate column
lm_coeff_df <- lm_coeff_df %>%
    filter(
        !grepl("const", variate)
    )
# if the log10p is inf then set to the max value
lm_coeff_df$log10p_value[is.infinite(lm_coeff_df$log10p_value)] <- max(
    lm_coeff_df$log10p_value[!is.infinite(lm_coeff_df$log10p_value)]
)

lm_coeff_df$Feature_type <- gsub(
    "RadialDistribution",
    "Radial\nDistibution",
    lm_coeff_df$Feature_type
)

# set the feature to feature + feature_number if scDINO in featurizer id
lm_coeff_df %>%
    mutate(
        feature = ifelse(
            grepl("scDINO", featurizer_id),
            paste0(feature, "_", feature_number),
            feature
        )
    ) -> lm_coeff_df


cp_and_scDINO_lm_coeff_df <- lm_coeff_df
cp_lm_coeff_df <- lm_coeff_df %>%
    filter(
        grepl("CP", featurizer_id)
    )
scDINO_lm_coeff_df <- lm_coeff_df %>%
    filter(
        grepl("scDINO", featurizer_id)
    )
list_of_dfs <- list(
    "both" = cp_and_scDINO_lm_coeff_df,
    "cp" = cp_lm_coeff_df,
    "scDINO" = scDINO_lm_coeff_df
)

for (df_name in names(list_of_dfs)) {
    df <- list_of_dfs[[df_name]]
    print(paste0("==============================="))
    print(paste0("Processing df: ", df_name))
    total_features <- length(unique(df$feature))
    print(paste0("Number of features: ", total_features))
    all_cell_count <- df %>%
        dplyr::filter(
            variate == "Cell count" & p_value_corrected < 0.05
        ) %>%
        dplyr::pull(feature)
    all_cell_count <- unique(all_cell_count)
    all_dose <- df %>%
        dplyr::filter(
            variate == "Dose" & p_value_corrected < 0.05
        ) %>%
        dplyr::pull(feature)
    all_dose <- unique(all_dose)
    all_time <- df %>%
        dplyr::filter(
            variate == "Time" & p_value_corrected < 0.05
        ) %>%
        dplyr::pull(feature)
    all_time <- unique(all_time)
    all_intersection <- Reduce(
        intersect,
        list(all_cell_count, all_dose, all_time)
    )
    # get the features that are unique to each variate
    cell_count_unique <- setdiff(
        all_cell_count,
        union(all_dose, all_time)
    )
    dose_unique <- setdiff(
        all_dose,
        union(all_cell_count, all_time)
    )
    time_unique <- setdiff(
        all_time,
        union(all_cell_count, all_dose)
    )

    cell_time_intersection <- setdiff(
        Reduce(
            intersect,
            list(all_cell_count, all_time)
        ),
        all_dose
    )
    cell_dose_intersection <- setdiff(
        Reduce(
            intersect,
            list(all_cell_count, all_dose)
        ),
        all_time
    )
    time_dose_intersection <- setdiff(
        Reduce(
            intersect,
            list(all_time, all_dose)
        ),
        all_cell_count
    )


    print(
        paste0(
            "Number of features significant for all three variates: ",
            length(all_intersection),
            " Percent of total: ",
            round(length(all_intersection)/total_features*100, 2), "%"
        )
    )
    print(
        paste0(
            "Number of features significant for all Time: ",
            length(all_time),
            " Percent of total: ",
            round(length(all_time)/total_features*100, 2), "%"
        )
    )
    print(
        paste0(
            "Number of features significant for all Dose: ",
            length(all_dose),
            " Percent of total: ",
            round(length(all_dose)/total_features*100, 2), "%"
        )
    )
    print(
        paste0(
            "Number of features significant for all Cell count: ",
            length(all_cell_count),
            " Percent of total: ",
            round(length(all_cell_count)/total_features*100, 2), "%"
        )
    )
    print(
        paste0(
            "Number of features significant for Cell count and Time: ",
            length(cell_time_intersection),
            " Percent of total: ",
            round(length(cell_time_intersection)/total_features*100, 2), "%"
        )
    )
    print(
        paste0(
            "Number of features significant for Cell count and Dose: ",
            length(cell_dose_intersection),
            " Percent of total: ",
            round(length(cell_dose_intersection)/total_features*100, 2), "%"
        )
    )
    print(
        paste0(
            "Number of features significant for Time and Dose: ",
            length(time_dose_intersection),
            " Percent of total: ",
            round(length(time_dose_intersection)/total_features*100, 2), "%"
        )
    )
    print(
        paste0(
            "Number of features significant for Cell count only: ",
            length(cell_count_unique),
            " Percent of total: ",
            round(length(cell_count_unique)/total_features*100, 2), "%"
        )
    )
    print(
        paste0(
            "Number of features significant for Dose only: ",
            length(dose_unique),
            " Percent of total: ",
            round(length(dose_unique)/total_features*100, 2), "%"
        )
    )
    print(
        paste0(
            "Number of features significant for Time only: ",
            length(time_unique),
            " Percent of total: ",
            round(length(time_unique)/total_features*100, 2), "%"
        )
    )
}

