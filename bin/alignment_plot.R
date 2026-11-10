library(ggplot2)
library(dplyr)
library(scales)
library(stringr)

theme_set(theme_bw())

args <- commandArgs(trailingOnly = TRUE)

tsv <- args[1]
workdir <- args[2]
prefix <- args[3]
datatype <- args[4]
skip_id_filter <- tolower(args[5])
total_reads_pre_filter <- as.integer(args[6])

setwd(workdir)

species_to_plot <- 14

df <- read.table(tsv, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

totalreads_int <- as.integer(sum(df$num_reads))

# adding commas for every thousanths place present in totalreads var
totalreads <- format(totalreads_int, big.mark = ",", scientific = FALSE)

# trim any leading spaces from 3rd column in dataframe
df$percent_reads <- sub("^\\s+", "", df$percent_reads)

# convert 3rd column to double
df$percent_reads <- as.double(df$percent_reads)

# grab the top 20 rows by total_species_count
unique_species <- length(unique(df$species))

if (unique_species > species_to_plot) {

	#	 top_species_names <- df %>%
	#		 group_by(species) %>%
	#		 summarise(total_reads = sum(num_reads, na.rm = TRUE)) %>%
	#		 slice_max(n = species_to_plot, order_by = total_reads) %>%
	#		 pull(species)
		
	#	 print(top_species_names)

	#	 top_species <- df %>%
	#		 filter(species %in% top_species_names)

	counts <- c()

	for (i in seq_len(nrow(df))) {

		species <- as.character(df$species[i])
		num_reads	<- df$num_reads[i]

		if (species %in% names(counts)) {

			counts[species] <- counts[species] + num_reads

		} else {

			counts[species] <- num_reads
			
		}
	}

	top_species_names <- names(head(sort(counts, decreasing = TRUE), species_to_plot))

	top_species <- df %>%
		filter(species %in% top_species_names)

	# adding special logic to include homo sapiens in the graph even if they're not in the top_species df - consequently kicks out the lowest count species in the df (so it won't be plotted)
	if ( 'homo sapiens' %in% tolower(top_species_names) ) {

		species_to_remove <- 'no_species_designated'

	} else {

		species_to_remove <- tail(top_species_names, 1)

	}

	# if we want to grab the top 20 species we want to make sure we also grab all ambiguity entries of each species, even if the read count for one ambiguity level is much lower than the other
	for (i in 1:nrow(df)) {

		row <- df[i, ]

		if ( (row$species %in% top_species$species | row$species == 'Homo sapiens') && row$species != species_to_remove) {

			if (exists("top_hits")) {

				top_hits <- rbind(top_hits, row)

			} else {

				top_hits <- row

			}

		}
	
	}

} else {

	top_hits <- df

}

top_hits_mod <- top_hits %>%
	arrange(species)

# adding spaces to really long species names to ensure proper string wrapping of axis labels during plotting stage
insert_space_every_25 <- function(s) {
	chars <- strsplit(s, "")[[1]]
	result <- character(0)
	pre_space_count <- 0
	count <- 0
	for (ch in chars) {
		if (ch == " ") {
			pre_space_count <- count
			count <- 0
			result <- c(result, ch)
		} else {
			if (count + pre_space_count == 23) {
				result <- c(result, " ")
				pre_space_count <- 0
				count <- 0
			}
			result <- c(result, ch)
			count <- count + 1
		}
	}
	paste(result, collapse = "")
}

top_hits_mod$species <- sapply(top_hits_mod$species, insert_space_every_25)

### This block will calculate the number and percent of "other" species classified by the "minimap2" algorithm and add them as a new row
total_percent <- sum(top_hits_mod$percent_reads)

if (total_percent < 99.9999 && unique_species > species_to_plot) {

	# reference original df to grab read count of species (not just top 20)
	read_count_all <- sum(df$num_reads)

	percent_other <- 100 - total_percent
	count_other <- read_count_all - sum(top_hits_mod$num_reads)

	new_row_other <- data.frame(
	species = "other",
		num_reads = count_other,
		percent_reads = percent_other,
		ambiguity = "none"
	)
	top_hits_mod <- rbind(top_hits_mod, new_row_other)

}

color_list <- c(
	"#FF7666" # Lighter Red
)

top_hits_mod$color <- NA

#setting up a separate dataframe that will make it easier to plot with ggplot
ggplot_df <- top_hits_mod

ggplot_df$color[ggplot_df$ambiguity == "none"] <- "#FF7666"
ggplot_df$color[ggplot_df$ambiguity == "single_genome"] <- "#ADD8E6"
ggplot_df$color[ggplot_df$ambiguity == "multi_genome"] <- "#D8BFD8"

ggplot_df$total_species_count <- NA
ggplot_df$total_species_percent <- NA
ggplot_df$offset <- NA

# creating samll df to see which species contains the highest read percentage regardless of ambiguity
total_count_by_species <- aggregate(num_reads ~ species, data = top_hits_mod, sum)
total_percent_by_species <- aggregate(percent_reads ~ species, data = top_hits_mod, sum)

# grabbing highest percent value to adjust data label positions in plot. The highest percent value will ultimately determine the scale of the y-axis, so our data label positions will be based on that
max_percent = max(total_percent_by_species$percent_reads, na.rm = TRUE)

# size factor is used to estimate the height of a single data label - the number 15 was determined via trial and error
size_factor = max_percent / 30

ggplot_df <- ggplot_df %>%
	arrange(species, factor(ambiguity, levels = c("none", "single_genome", "multi_genome")))

# this for loop will determine the y-axis position of each data label to ensure there is no overlap
for (i in 1:nrow(ggplot_df)) {

	row <- ggplot_df[i, ]

	ggplot_df[i, "total_species_count"] <- total_count_by_species$num_reads[total_count_by_species$species == row$species]

	ggplot_df[i, "total_species_percent"] <- total_percent_by_species$percent_reads[total_percent_by_species$species == row$species]

	if (row$ambiguity == "none") {

		ggplot_df[i, "offset"] <- (total_percent_by_species$percent_reads[total_percent_by_species$species == row$species] - ( total_percent_by_species$percent_reads[total_percent_by_species$species == row$species] - row$percent_reads)) / 2

	} else if (row$ambiguity == "single_genome") {

		none_ambiguity_percent <- ggplot_df$percent_reads[ggplot_df$species == row$species & ggplot_df$ambiguity == "none"]
		if (length(none_ambiguity_percent) == 0) none_ambiguity_percent <- 0

		# center of the single_genome segment
		vertical_pos <- none_ambiguity_percent + row$percent_reads / 2

		# if there is a none segment and its label would be too close, nudge up
		if (none_ambiguity_percent > 0 && vertical_pos - none_ambiguity_percent / 2 < size_factor) {
			vertical_pos <- none_ambiguity_percent / 2 + size_factor
		}

		ggplot_df[i, "offset"] <- vertical_pos

	} else {

		none_ambiguity_percent = ggplot_df$percent_reads[ggplot_df$species == row$species & ggplot_df$ambiguity == "none"]
		single_ambiguity_percent = ggplot_df$percent_reads[ggplot_df$species == row$species & ggplot_df$ambiguity == "single_genome"]
		single_ambiguity_offset = ggplot_df$offset[ggplot_df$species == row$species & ggplot_df$ambiguity == "single_genome"]

		vertical_pos = total_percent_by_species$percent_reads[total_percent_by_species$species == row$species] - (row$percent_reads / 2)

		# determines if 'single_genome' and 'none' ambiguity classifiers exists for this species; this is required to prevent a 'length of zero' error in the subsequent if statements
		if (length(single_ambiguity_offset) == 0 & length(none_ambiguity_percent) == 0) {

			vertical_pos = vertical_pos

		} else if (length(single_ambiguity_offset) == 0 & length(none_ambiguity_percent) != 0) {

			if (none_ambiguity_percent >= size_factor) {

				vertical_pos = vertical_pos
			
			} else {

				vertical_pos = none_ambiguity_percent + size_factor

			}

		} else if (length(single_ambiguity_offset) != 0 & length(none_ambiguity_percent) == 0) {

				if (single_ambiguity_offset + size_factor <= vertical_pos) {

					vertical_pos = vertical_pos
				
				} else {

					vertical_pos = single_ambiguity_offset + size_factor

				}

		} else if (single_ambiguity_offset + size_factor <= vertical_pos) {

			vertical_pos = vertical_pos
		
		} else {

			vertical_pos = single_ambiguity_offset + size_factor

		}

		ggplot_df[i, "offset"] <- vertical_pos

	}

}

# format subtitle depending on whether the read set was filtered by tax ID before hand
if (skip_id_filter == 'false') {

	total_reads_pre_filter_as_text <- format(total_reads_pre_filter, big.mark = ",", scientific = FALSE)

	post_filter_percent <- round((totalreads_int / total_reads_pre_filter) * 100, 2)

	subtitle_text <- paste0("Total ", datatype, " (Pre Tax ID Filter): ", total_reads_pre_filter_as_text, "	 ", "Total ", datatype, " (Post Tax ID Filter): ", totalreads, " (", post_filter_percent, "%)")

} else {

subtitle_text <- paste0("Total ", datatype, ": ", totalreads)

}

# ordering species alphabetically while ignoring capitalization of the first letter and setting levels to ambiguity vector
ggplot_df <- ggplot_df %>%
	arrange(tolower(species), factor(ambiguity, levels = c("none", "single_genome", "multi_genome")))

# ggplot_df$ambiguity <- factor(ggplot_df$ambiguity, levels = rev(c("none", "single_genome", "multi_genome")))

# separate out 'unclassified' and 'other' rows via a for loop, then reassemble with them at the end
for (i in 1:nrow(ggplot_df)) {

	row <- ggplot_df[i, ]

	if (row$species == "unclassified" | row$species == "other" | tolower(row$species) == "homo sapiens") {

		if (exists("special_rows")) {

			special_rows <- rbind(special_rows, row)

		} else {

			special_rows <- row

		}

	} else {

		if (exists("regular_rows")) {

			regular_rows <- rbind(regular_rows, row)

		} else {

			regular_rows <- row

		}

	}

}

# if special rows exist, order them as unclassified then other, and reattach at the end
if (exists("special_rows")) {

	special_rows <- special_rows %>%
		arrange(factor(species, levels = c("Homo sapiens", "homo sapiens", "unclassified", "other")))

	ggplot_df <- rbind(regular_rows, special_rows)

} else {

	ggplot_df <- regular_rows

}

# force ggplot2 to use the order of species as they appear in ggplot_df by converting species column to a factor with levels set in the desired order
# scale_x_discrete is a little wonky looking, but it will wrap x-axis labels at 25 chars and cut them off at 50 chars - prevents the graph from being scaled down to fit the name of the organism into the png 
ggplot_df$species <- factor(ggplot_df$species, levels = unique(ggplot_df$species))
ggplot_df$ambiguity <- factor(ggplot_df$ambiguity, levels = c("multi_genome", "single_genome", "none"))

plot_1 <- ggplot(ggplot_df, aes(x = species, group = ambiguity, y = percent_reads)) +
	geom_col(aes(fill = ambiguity), color = "black") +
	theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10), plot.margin = margin(t = 10, r = 5, b = 5, l = 15)) +
	scale_x_discrete(labels = function(x) label_wrap(25)(str_trunc(x, width = 50, ellipsis = "..."))) +
	scale_fill_manual(values = setNames(ggplot_df$color, ggplot_df$ambiguity)) +
	labs(
		title = paste(prefix, sub("s$", "",datatype), "Taxonomic Distribution"),
		subtitle = subtitle_text,
		x = "",
		y = paste0("Percent Mapped")
	) +
	geom_label(
		data = ggplot_df,
		size = 2.8,
		aes(label=paste0(sprintf("%.2f", percent_reads), "%", " | ", num_reads), y = offset, fill = factor(ambiguity)),
		label.padding = unit(0.15, "lines"),
		show.legend	= FALSE
	)

ggsave("plot.png", width = 16, height = 7, dpi = 150)
