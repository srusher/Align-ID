################################################################################
# Taxonomic distribution plot
#
# Produces a stacked bar chart (plot.png) of the percent of reads/read pairs
# mapped to each species, split by ambiguity class (none / single_genome /
# multi_genome). Each bar has:
#   - data labels for every ambiguity segment ("percent | read count")
#   - a small donut chart + text above it showing genome coverage (%)
#
# Usage:
#   Rscript plot.R <tsv> <workdir> <prefix> <datatype> <skip_id_filter> \
#                  <total_reads_pre_filter> <coverage_tsv>
#
# Args:
#   tsv                    - read classification table (species, num_reads,
#                            percent_reads, ambiguity)
#   workdir                - directory to work in / write plot.png to
#   prefix                 - sample name used in the plot title
#   datatype               - e.g. "Reads" or "Read Pairs" (used in title/subtitle)
#   skip_id_filter         - "true"/"false"; whether a tax ID filter was skipped
#   total_reads_pre_filter - total reads before tax ID filtering (for subtitle)
#   coverage_tsv           - per-species coverage table (needs 'species' and
#                            'genome_cov(%)' columns)
################################################################################

library(ggplot2)
library(dplyr)
library(scales)   # label_wrap() for x-axis labels
library(stringr)  # str_trunc() for x-axis labels
library(grid)     # low-level grobs used to draw the coverage donuts

theme_set(theme_bw())

# ---- command line arguments --------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

tsv <- args[1]
workdir <- args[2]
prefix <- args[3]
datatype <- args[4]
skip_id_filter <- tolower(args[5])
total_reads_pre_filter <- as.integer(args[6])
coverage_tsv <- args[7]

setwd(workdir)

# max number of species (bars) to show before collapsing the rest into "other"
species_to_plot <- 14

# ---- load classification data ------------------------------------------------
df <- read.table(tsv, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# total reads across all species/ambiguity rows (used in the subtitle)
totalreads_int <- as.integer(sum(df$num_reads))

# adding commas for every thousanths place present in totalreads var
totalreads <- format(totalreads_int, big.mark = ",", scientific = FALSE)

# trim any leading spaces from 3rd column in dataframe
df$percent_reads <- sub("^\\s+", "", df$percent_reads)

# convert 3rd column to double
df$percent_reads <- as.double(df$percent_reads)

# number of distinct species in the input; decides whether we need to filter to the top N
unique_species <- length(unique(df$species))

# ---- select the species to plot ----------------------------------------------
if (unique_species > species_to_plot) {

	# (old dplyr-based approach to finding the top species, replaced by the loop below)
	#	 top_species_names <- df %>%
	#		 group_by(species) %>%
	#		 summarise(total_reads = sum(num_reads, na.rm = TRUE)) %>%
	#		 slice_max(n = species_to_plot, order_by = total_reads) %>%
	#		 pull(species)

	#	 print(top_species_names)

	#	 top_species <- df %>%
	#		 filter(species %in% top_species_names)

	# named vector: total read count per species, summed across ambiguity classes
	counts <- c()

	for (i in seq_len(nrow(df))) {

		species <- as.character(df$species[i])
		num_reads	<- df$num_reads[i]

		if (species %in% names(counts)) {

			# species already seen - add this row's reads to its running total
			counts[species] <- counts[species] + num_reads

		} else {

			# first time seeing this species
			counts[species] <- num_reads

		}
	}

	# names of the top N species by total read count
	top_species_names <- names(head(sort(counts, decreasing = TRUE), species_to_plot))

	# all rows (every ambiguity class) belonging to those top species
	top_species <- df %>%
		filter(species %in% top_species_names)

	# adding special logic to include homo sapiens in the graph even if they're not in the top_species df - consequently kicks out the lowest count species in the df (so it won't be plotted)
	if ( 'homo sapiens' %in% tolower(top_species_names) ) {

		# human is already in the top N, so nothing needs to be kicked out
		species_to_remove <- 'no_species_designated'

	} else {

		# human is forced in below, so drop the lowest-ranked species to keep the bar count constant
		species_to_remove <- tail(top_species_names, 1)

	}

	# if we want to grab the top N species we want to make sure we also grab all ambiguity entries of each species, even if the read count for one ambiguity level is much lower than the other
	for (i in 1:nrow(df)) {

		row <- df[i, ]

		# keep rows for top species (plus Homo sapiens), minus the species being kicked out
		if ( (row$species %in% top_species$species | row$species == 'Homo sapiens') && row$species != species_to_remove) {

			if (exists("top_hits")) {

				top_hits <- rbind(top_hits, row)

			} else {

				top_hits <- row

			}

		}

	}

} else {

	# few enough species that everything gets plotted
	top_hits <- df

}

top_hits_mod <- top_hits %>%
	arrange(species)

# ---- axis label formatting ---------------------------------------------------
# adding spaces to really long species names to ensure proper string wrapping of axis labels during plotting stage
# (walks the name character by character and inserts a space once ~23 chars have gone by
# without one; the counter tracks the length of the previous word so wrapping stays balanced)
insert_space_every_25 <- function(s) {
	chars <- strsplit(s, "")[[1]]
	result <- character(0)
	pre_space_count <- 0   # length of the previous word
	count <- 0             # chars since the last space
	for (ch in chars) {
		if (ch == " ") {
			# natural space: remember how long the word was and reset
			pre_space_count <- count
			count <- 0
			result <- c(result, ch)
		} else {
			# word is running too long - force a break here
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

# if the plotted species don't account for 100% of reads, lump the remainder into "other"
if (total_percent < 99.9999 && unique_species > species_to_plot) {

	# reference original df to grab read count of species (not just top N)
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

# ---- genome coverage lookup --------------------------------------------------
# read the coverage table; check.names = FALSE keeps the 'genome_cov(%)' column name intact
cov_df <- read.delim(coverage_tsv, sep = "\t", header = TRUE,
                     stringsAsFactors = FALSE, check.names = FALSE,
                     quote = "", comment.char = "")

# keep only species + numeric coverage. Species names get the same space-insertion
# as above so they match top_hits_mod$species; the "%" sign is stripped from coverage values
cov_df <- data.frame(
	species = sapply(cov_df$species, insert_space_every_25),
	genome_cov = as.numeric(sub("%", "", cov_df[["genome_cov(%)"]])),
	stringsAsFactors = FALSE
) %>%
	distinct(species, .keep_all = TRUE)

# text formatter for coverage values (NA -> blank, tiny non-zero -> "<0.01% cov")
# NOTE: currently unused for the donut text (which uses sprintf directly), but used for cov_label below
format_cov <- function(x) {
	ifelse(is.na(x), "",
		ifelse(x > 0 & x < 0.01, "<0.01% cov",
			paste0(sprintf("%.2f", x), "% cov")))
}

# attach coverage to every row; "other" and unclassified rows get NA
top_hits_mod <- top_hits_mod %>%
	left_join(cov_df, by = "species")

# ---- fill colors per ambiguity class -----------------------------------------
color_list <- c(
	"#FF7666" # Lighter Red
)

top_hits_mod$color <- NA

#setting up a separate dataframe that will make it easier to plot with ggplot
ggplot_df <- top_hits_mod

ggplot_df$color[ggplot_df$ambiguity == "none"] <- "#FF7666"          # red
ggplot_df$color[ggplot_df$ambiguity == "single_genome"] <- "#ADD8E6" # light blue
ggplot_df$color[ggplot_df$ambiguity == "multi_genome"] <- "#D8BFD8"  # light purple

# placeholder columns filled in by the label-positioning loop below
ggplot_df$total_species_count <- NA
ggplot_df$total_species_percent <- NA
ggplot_df$offset <- NA   # y position of each segment's data label

# creating samll df to see which species contains the highest read percentage regardless of ambiguity
total_count_by_species <- aggregate(num_reads ~ species, data = top_hits_mod, sum)
total_percent_by_species <- aggregate(percent_reads ~ species, data = top_hits_mod, sum)

# grabbing highest percent value to adjust data label positions in plot. The highest percent value will ultimately determine the scale of the y-axis, so our data label positions will be based on that
max_percent = max(total_percent_by_species$percent_reads, na.rm = TRUE)

# size factor is used to estimate the height of a single data label (in y-axis units) - the divisor of 30 was determined via trial and error
size_factor = max_percent / 30

# order rows within each species: none -> single_genome -> multi_genome (bottom to top of the stack)
ggplot_df <- ggplot_df %>%
	arrange(species, factor(ambiguity, levels = c("none", "single_genome", "multi_genome")))

# this for loop will determine the y-axis position of each data label to ensure there is no overlap
# (rows are processed bottom-to-top within a species, so each label can look at the
# offsets of the segments below it and nudge itself upward if it would collide)
for (i in 1:nrow(ggplot_df)) {

	row <- ggplot_df[i, ]

	# record per-species totals on each row
	ggplot_df[i, "total_species_count"] <- total_count_by_species$num_reads[total_count_by_species$species == row$species]

	ggplot_df[i, "total_species_percent"] <- total_percent_by_species$percent_reads[total_percent_by_species$species == row$species]

	if (row$ambiguity == "none") {

		# bottom segment: label goes in the middle of the segment
		ggplot_df[i, "offset"] <- (total_percent_by_species$percent_reads[total_percent_by_species$species == row$species] - ( total_percent_by_species$percent_reads[total_percent_by_species$species == row$species] - row$percent_reads)) / 2

	} else if (row$ambiguity == "single_genome") {

		# height of the 'none' segment beneath this one (0 if the species has none)
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

		# multi_genome (top segment): look up the heights/label positions of the segments below it
		none_ambiguity_percent = ggplot_df$percent_reads[ggplot_df$species == row$species & ggplot_df$ambiguity == "none"]
		single_ambiguity_percent = ggplot_df$percent_reads[ggplot_df$species == row$species & ggplot_df$ambiguity == "single_genome"]
		single_ambiguity_offset = ggplot_df$offset[ggplot_df$species == row$species & ggplot_df$ambiguity == "single_genome"]

		# default position: center of the top segment (top of stack minus half its height)
		vertical_pos = total_percent_by_species$percent_reads[total_percent_by_species$species == row$species] - (row$percent_reads / 2)

		# determines if 'single_genome' and 'none' ambiguity classifiers exists for this species; this is required to prevent a 'length of zero' error in the subsequent if statements
		if (length(single_ambiguity_offset) == 0 & length(none_ambiguity_percent) == 0) {

			# only a multi_genome segment exists - nothing to avoid
			vertical_pos = vertical_pos

		} else if (length(single_ambiguity_offset) == 0 & length(none_ambiguity_percent) != 0) {

			# 'none' + 'multi' only: nudge up if the 'none' segment is too short to leave room
			if (none_ambiguity_percent >= size_factor) {

				vertical_pos = vertical_pos

			} else {

				vertical_pos = none_ambiguity_percent + size_factor

			}

		} else if (length(single_ambiguity_offset) != 0 & length(none_ambiguity_percent) == 0) {

				# 'single' + 'multi' only: keep one label-height above the single label
				if (single_ambiguity_offset + size_factor <= vertical_pos) {

					vertical_pos = vertical_pos

				} else {

					vertical_pos = single_ambiguity_offset + size_factor

				}

		} else if (single_ambiguity_offset + size_factor <= vertical_pos) {

			# all three classes exist and there's already enough room
			vertical_pos = vertical_pos

		} else {

			# all three classes exist but labels would overlap - push up above the single label
			vertical_pos = single_ambiguity_offset + size_factor

		}

		ggplot_df[i, "offset"] <- vertical_pos

	}

}

# ---- subtitle ----------------------------------------------------------------
# format subtitle depending on whether the read set was filtered by tax ID before hand
if (skip_id_filter == 'false') {

	total_reads_pre_filter_as_text <- format(total_reads_pre_filter, big.mark = ",", scientific = FALSE)

	# what fraction of the original reads survived the tax ID filter
	post_filter_percent <- round((totalreads_int / total_reads_pre_filter) * 100, 2)

	subtitle_text <- paste0("Total ", datatype, " (Pre Tax ID Filter): ", total_reads_pre_filter_as_text, "	 ", "Total ", datatype, " (Post Tax ID Filter): ", totalreads, " (", post_filter_percent, "%)")

} else {

subtitle_text <- paste0("Total ", datatype, ": ", totalreads)

}

# ---- x-axis ordering ---------------------------------------------------------
# ordering species alphabetically while ignoring capitalization of the first letter and setting levels to ambiguity vector
ggplot_df <- ggplot_df %>%
	arrange(tolower(species), factor(ambiguity, levels = c("none", "single_genome", "multi_genome")))

# ggplot_df$ambiguity <- factor(ggplot_df$ambiguity, levels = rev(c("none", "single_genome", "multi_genome")))

# separate out 'unclassified' and 'other' rows via a for loop, then reassemble with them at the end
# (Homo sapiens is also treated as a "special" row so it sits at the far right of the plot)
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

# if special rows exist, order them as Homo sapiens, then unclassified, then other, and reattach at the end
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
# stacking order / legend order: multi_genome on top, none on the bottom
ggplot_df$ambiguity <- factor(ggplot_df$ambiguity, levels = c("multi_genome", "single_genome", "none"))

# ---- per-species data for the coverage donuts --------------------------------
# one row per species: bar height, where the highest data label ends, and coverage
cov_label_df <- ggplot_df %>%
	group_by(species) %>%
	summarise(
		total_percent = sum(percent_reads),                     # full height of the stacked bar
		label_top     = max(offset) + size_factor * 0.75,   # top edge of highest data label
		genome_cov    = first(genome_cov),
		.groups = "drop"
	) %>%
	mutate(
		# donut sits above whichever is higher: the bar or its tallest data label (avoids overlap on tiny bars)
		anchor_y  = pmax(total_percent, label_top),
		cov_label = format_cov(genome_cov)
	)

# builds a small donut chart (gray ring + colored arc proportional to coverage) with a
# text label above it. Everything is sized in inches and anchored to the bottom of its
# drawing region, so the donut stays perfectly circular regardless of the plot's axis scales
make_donut <- function(cov, label,
                       r_out = 0.16,   # outer radius, inches
                       r_in  = 0.09,   # inner radius (hole), inches
                       pad   = 0.04,   # gap between bar top and donut, inches
                       fill_col = "#019114",
                       bg_col   = "#D9D9D9") {

	# ring segment from angle t0 to t1 (radians, clockwise from 12 o'clock)
	# drawn as a polygon: outer arc forward, then inner arc backward
	band <- function(t0, t1, fill) {
		n <- max(2, ceiling(100 * (t1 - t0) / (2 * pi)))   # ~100 points per full circle
		t <- seq(t0, t1, length.out = n)
		x <- c(r_out * sin(t), rev(r_in * sin(t)))
		y <- c(r_out * cos(t), rev(r_in * cos(t)))
		polygonGrob(
			x = unit(0.5, "npc") + unit(x, "in"),
			y = unit(0, "npc") + unit(pad + r_out + y, "in"),
			gp = gpar(fill = fill, col = NA)
		)
	}

	# coverage as a 0-1 fraction (clamped to 0-100%)
	frac <- min(max(cov, 0), 100) / 100

	# full gray ring as the background track
	kids <- gList(band(0, 2 * pi, bg_col))
	# colored arc on top, covering the covered fraction of the ring
	if (frac > 0) kids <- gList(kids, band(0, 2 * pi * frac, fill_col))

	# coverage text centered just above the donut
	kids <- gList(
		kids,
		textGrob(
			label,
			x = unit(0.5, "npc"),
			y = unit(0, "npc") + unit(pad + 2 * r_out + 0.03, "in"),
			vjust = 0,
			gp = gpar(fontsize = 8, fontface = "bold")
		)
	)

	gTree(children = kids)
}

# one annotation_custom layer per species: places that species' donut in the x-slot
# of its bar (x = factor position +/- 0.5), starting at its anchor height
donut_layers <- lapply(seq_len(nrow(cov_label_df)), function(i) {

	r <- cov_label_df[i, ]
	if (is.na(r$genome_cov)) return(NULL)   # no donut for other/unclassified

	annotation_custom(
		grob = make_donut(r$genome_cov, sprintf("%.2f%%", r$genome_cov)),
		xmin = as.integer(r$species) - 0.5,
		xmax = as.integer(r$species) + 0.5,
		ymin = r$anchor_y,
		ymax = r$anchor_y + size_factor
	)
})

# legend key: a donut at ~75% fill, sized in inches so it fits the standard key box
draw_key_donut <- function(data, params, size,
                           r_out = 0.10, r_in = 0.055, frac = 0.75,
                           fill_col = "#019114", bg_col = "#D9D9D9") {

	# ring segment from t0 to t1 (radians, clockwise from 12 o'clock), centered in the key
	band <- function(t0, t1, fill) {
		n <- max(2, ceiling(100 * (t1 - t0) / (2 * pi)))
		t <- seq(t0, t1, length.out = n)
		x <- c(r_out * sin(t), rev(r_in * sin(t)))
		y <- c(r_out * cos(t), rev(r_in * cos(t)))
		polygonGrob(
			x = unit(0.5, "npc") + unit(x, "in"),
			y = unit(0.5, "npc") + unit(y, "in"),
			gp = gpar(fill = fill, col = NA)
		)
	}

	gTree(children = gList(band(0, 2 * pi, bg_col), band(0, 2 * pi * frac, fill_col)))
}

# drop the NULLs (species without coverage) so ggplot can add the list of layers
donut_layers <- Filter(Negate(is.null), donut_layers)

# ---- build the plot ----------------------------------------------------------
plot_1 <- ggplot(ggplot_df, aes(x = species, group = ambiguity, y = percent_reads)) +
	# stacked bars, filled by ambiguity class
	geom_col(aes(fill = ambiguity), color = "black") +
	# angled x-axis labels; margins leave room for the y-axis title
	theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10), plot.margin = margin(t = 10, r = 5, b = 5, l = 15)) +
	# wrap species names at 25 chars, truncate at 50
	scale_x_discrete(labels = function(x) label_wrap(25)(str_trunc(x, width = 50, ellipsis = "..."))) +
	# map each ambiguity level to its color
	scale_fill_manual(values = setNames(ggplot_df$color, ggplot_df$ambiguity)) +
	labs(
		title = paste(prefix, sub("s$", "",datatype), "Taxonomic Distribution"),
		subtitle = subtitle_text,
		x = "",
		y = paste0("Percent Mapped")
	) +
	# data labels ("percent | read count") at the offsets computed earlier
	geom_label(
		data = ggplot_df,
		size = 2.8,
		aes(label=paste0(sprintf("%.2f", percent_reads), "%", " | ", num_reads), y = offset, fill = factor(ambiguity)),
		label.padding = unit(0.15, "lines"),
		show.legend	= FALSE
	) +
	# coverage donuts above each bar
	donut_layers +
	geom_point(
		data = data.frame(species = ggplot_df$species[1], percent_reads = NA_real_, cov = "Genome coverage"),
		aes(x = species, y = percent_reads, color = cov),
		inherit.aes = FALSE,
		key_glyph = draw_key_donut,
		na.rm = TRUE,
		show.legend = c(color = TRUE, fill = FALSE)
	) +
	scale_color_manual(name = NULL, values = c("Genome coverage" = "#019114")) +
	scale_y_continuous(expand = expansion(mult = c(0.02, 0.12)))

# NOTE: ggsave() with no plot argument saves the last plot displayed (plot_1 here)
ggsave("plot.png", width = 16, height = 7, dpi = 150)