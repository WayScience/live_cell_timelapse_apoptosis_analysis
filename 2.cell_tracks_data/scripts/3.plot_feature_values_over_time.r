suppressPackageStartupMessages(suppressWarnings({
    library("ggplot2")
    library(dplyr)
    library(tidyr)
    library(ComplexHeatmap)
    library(tibble)
    library(RColorBrewer)
    library(scales)
    library(circlize)
}))


profile_file_path <- file.path("../../data/CP_feature_select/profiles/features_selected_profile.parquet")
figure_path <- file.path("../figures/")
if (!dir.exists(figure_path)) {
    dir.create(figure_path, recursive = TRUE)
}

profile <- arrow::read_parquet(profile_file_path, col_select = c("Metadata_dose", "Metadata_Well", "Metadata_FOV"))

# arrange by well_fov
profile <- profile %>% arrange(Metadata_Well) %>%
    select(Metadata_Well, Metadata_dose)
# remove duplicates
profile <- profile %>%
    group_by(Metadata_Well) %>%
    summarise(Metadata_dose = unique(Metadata_dose)) %>%
    ungroup()

profile_file_path <- file.path("../../data/CP_aggregated/profiles/aggregated_profile.parquet")
df <- arrow::read_parquet(profile_file_path) %>% arrange(Metadata_Well)
df <- df %>%
    left_join(profile, by = c("Metadata_Well" = "Metadata_Well"), relationship ="many-to-many")
# transform the data to standard scalar (-1, 1) format
for (i in 1:ncol(df)) {
    # make sure the column is not metadata
    if (grepl("Metadata_", colnames(df)[i])) {
        next
    }
    if (is.numeric(df[[i]])) {
        df[[i]] <- rescale(df[[i]], to = c(-1, 1))
    }
}
# map each of the Time points to the actual timepoint
df$Metadata_Time <- as.numeric(df$Metadata_Time * 60 /2 )
head(df)

# complex heatmap does not compare across heatmaps the scale so we must set it manually
# for more information see:
# https://github.com/jokergoo/EnrichedHeatmap/issues/7

# we will set the color scale the same way that ComplexHeatmap does
# The automatically generated colors map from the minus and plus 99^th of
# the absolute values in the matrix.


global_across_dose_99th_min <- df %>%
    select(-Metadata_Well, -Metadata_dose, -Metadata_Time) %>%
    summarise(across(everything(), ~ quantile((.), 0.01, na.rm = TRUE))) %>%
    unlist() %>%
    min(na.rm = TRUE)
global_across_dose_99th_max <- df %>%
    select(-Metadata_Well, -Metadata_dose, -Metadata_Time) %>%
    summarise(across(everything(), ~ quantile((.), 0.99, na.rm = TRUE))) %>%
    unlist() %>%
    max(na.rm = TRUE)

print(global_across_dose_99th_min)
print(global_across_dose_99th_max)
col_fun = circlize::colorRamp2(c(global_across_dose_99th_min, 0, global_across_dose_99th_max), c("blue","white", "red"))

# get the list of features
features <- colnames(df)
features <- features[!features %in% c("Metadata_Well", "Metadata_dose", "Metadata_Time")]
features <- as.data.frame(features)
# split the features by _ into multiple columns
features <- features %>%
    separate(features, into = c("Compartment", "Measurement", "Metric", "Extra", "Extra1", "Extra2", "Extra3"), sep = "_", extra = "merge", fill = "right")
# clean up the features columns
# if Extra is NA then replace with None
features$Extra[is.na(features$Extra)] <- "None"
# if extra is a number then replace with None
features$Extra[grepl("^[0-9]+$", features$Extra)] <- "None"
# replace all other NAs with None
features$Extra1[is.na(features$Extra1)] <- "None"
features$Extra2[is.na(features$Extra2)] <- "None"
# change extra to None if X or Y
features$Extra[features$Extra == "X"] <- "None"
features$Extra[features$Extra == "Y"] <- "None"
# drop the Adjacent channel
features$Extra[features$Extra == "Adjacent"] <- "None"
# if extra1 is 488 then add extra2 to Extra1
features$Extra1[features$Extra1 == "488"] <- paste0(features$Extra1[features$Extra1 == "488"], "_", features$Extra2[features$Extra1 == "488"])
# if extra1 id CL then add extra1 to Extra
features$Extra[features$Extra == "CL"] <- paste0(features$Extra[features$Extra == "CL"], "_", features$Extra1[features$Extra == "CL"])

features <- features %>%
    rename(Channel = Extra) %>%
    select(-Extra1, -Extra2)
# rename channel names to replace "_" with " "
features$Channel <- gsub("CL_488_1", "CL 488_1", features$Channel)
features$Channel <- gsub("CL_488_2", "CL 488_2", features$Channel)
features$Channel <- gsub("CL_561", "CL 561", features$Channel)

# time color function
time_col_fun = colorRamp2(
    c(min(unique(df$Metadata_Time)), max(unique(df$Metadata_Time))), c("white", "purple")
    )

column_anno <- HeatmapAnnotation(
    Time = unique(df$Metadata_Time),
    show_legend = TRUE,
    annotation_name_gp = gpar(fontsize = 2),
    annotation_legend_param = list(
        title_position = "topcenter",
        title_gp = gpar(fontsize = 16, angle = 0, fontface = "bold", hjust = 1.0),
        labels_gp = gpar(fontsize = 16,
        title = gpar(fontsize = 16))),
    col = list(
        Time = time_col_fun
    )
)

# compartment row annotation
row_compartment = rowAnnotation(
    Object = features$Compartment,
        show_legend = TRUE,
    # change the legend titles
    annotation_legend_param = list(
        title_position = "topcenter",
        title_gp = gpar(fontsize = 16, angle = 0, fontface = "bold", hjust = 1.0),
        labels_gp = gpar(fontsize = 16,
        title = gpar(fontsize = 16))),
    annotation_name_side = "bottom",
    annotation_name_gp = gpar(fontsize = 16),
    # color
    col = list(
        Object = c(
            "Cells" = "#B000B0",
            "Cytoplasm" = "#00D55B",
            "Nuclei" = "#0000AB"
            )
    )
)
row_measurement = rowAnnotation(
    FeatureType = features$Measurement,
           annotation_legend_param = list(
        title_position = "topcenter",
        title_gp = gpar(fontsize = 16, angle = 0, fontface = "bold", hjust = 0.5),
        labels_gp = gpar(fontsize = 16,
        title = gpar(fontsize = 16))),
    annotation_name_side = "bottom",
    annotation_name_gp = gpar(fontsize = 16),
    col = list(
            FeatureType = c(
            "AreaShape" = brewer.pal(8, "Paired")[1],
            "Correlation" = brewer.pal(8, "Paired")[2],
            "Granularity" = brewer.pal(8, "Paired")[3],
            "Intensity" = brewer.pal(8, "Paired")[4],
            "Location" = brewer.pal(8, "Paired")[5],
            "Neighbors" =  brewer.pal(8, "Paired")[6],
            "RadialDistribution" = brewer.pal(8, "Paired")[7],
            "Texture" = brewer.pal(8, "Paired")[8]
        )
    ),
    show_legend = TRUE
)
row_channel = rowAnnotation(
    Channel = features$Channel,
        annotation_legend_param = list(
        title_position = "topcenter",
        title_gp = gpar(fontsize = 16, angle = 0, fontface = "bold", hjust = 0.5),
        labels_gp = gpar(fontsize = 16,
        # make annotation bar text bigger
        legend = gpar(fontsize = 16),
        annotation_name = gpar(fontsize = 16),
        legend_height = unit(20, "cm"),
        legend_width = unit(1, "cm"),
        # make legend taller
        legend_height = unit(10, "cm"),
        legend_width = unit(1, "cm"),
        legend_key = gpar(fontsize = 16)
        )
    ),



    annotation_name_side = "bottom",
    # make font size bigger
    annotation_name_gp = gpar(fontsize = 16),
    col = list(
    Channel = c(
            "DNA" = "#0000AB",
            "CL 488_1" = "#B000B0",
            "CL 488_2" = "#00D55B",
            "CL 561" = "#FFFF00",
            "None" = "#B09FB0")
    )
)
row_annotations = c(row_compartment, row_measurement, row_channel)

list_of_mats_for_heatmaps <- list()
list_of_heatmaps <- list()
heatmap_list <- NULL
ht_opt(RESET = TRUE)
df$Metadata_dose <- as.numeric(df$Metadata_dose)
for (dose in unique(df$Metadata_dose)) {
    # check if the last in the number of doses

    # get the first dose
    single_dose_df <- df %>%
        filter(Metadata_dose == dose) %>%
        group_by(Metadata_Time) %>%
        select(-Metadata_Well, -Metadata_dose) %>%
        summarise(across(everything(), ~ mean(., na.rm = TRUE))) %>%
        ungroup()

    # sort the columns by Metadata_Time
    single_dose_df <- single_dose_df %>%
        select(Metadata_Time, everything()) %>%
        arrange(Metadata_Time)

    mat <- t(as.matrix(single_dose_df))

    colnames(mat) <- single_dose_df$Metadata_Time
    mat <- mat[-1,]

    if (dose == max(unique(df$Metadata_dose))) {

        heatmap_plot <- Heatmap(
            mat,
            col = col_fun,
            show_row_names = FALSE,
            show_column_names = FALSE,
            cluster_columns = FALSE,
            column_names_gp = gpar(fontsize = 16), # Column name label formatting
            row_names_gp = gpar(fontsize = 14),

            show_heatmap_legend = TRUE,
            heatmap_legend_param = list(
                        title = "Feature\nValue",
                        title_position = "topcenter",
                        # direction = "horizontal",
                        title_gp = gpar(fontsize = 16, angle = 0, fontface = "bold", hjust = 1.0),
                        labels_gp = gpar(fontsize = 16),
                        legend_height = unit(4, "cm"),
                        legend_width = unit(3, "cm"),
                        annotation_legend_side = "bottom"
                        ),
            row_dend_width = unit(2, "cm"),
            column_title = paste0("Dose: ", dose," uM"),
            # add the row annotations
            right_annotation = row_annotations,
            top_annotation = column_anno
        )
    } else {
        heatmap_plot <- Heatmap(
            mat,
            col = col_fun,
            show_row_names = FALSE,
            cluster_columns = FALSE,
            show_column_names = FALSE,

            column_names_gp = gpar(fontsize = 16), # Column name label formatting
            row_names_gp = gpar(fontsize = 14),

            show_heatmap_legend = FALSE,
            heatmap_legend_param = list(
                        title = "Feature\nValue",
                        title_position = "topcenter",
                        title_gp = gpar(fontsize = 16, angle = 0, fontface = "bold", hjust = 1.0),
                        labels_gp = gpar(fontsize = 16),
                        legend_height = unit(4, "cm"),
                        legend_width = unit(3, "cm"),
                        annotation_legend_side = "bottom"
                        ),
            row_dend_width = unit(2, "cm"),
            column_title = paste0("Dose: ", dose," uM"),
            top_annotation = column_anno,
        )
    }
    # add the heatmap to the list
    heatmap_list <- heatmap_list + heatmap_plot
}

width <- 20
height <- 12
options(repr.plot.width = width, repr.plot.height = height)
ht_opt$message = FALSE
png(filename = paste0(figure_path, "filtered_features.png"), width = width, height = height, units = "in", res = 600)
draw(
    heatmap_list,
    merge_legends = TRUE,

    heatmap_legend_side = "right",
    annotation_legend_side = "bottom"
)
dev.off()
