install.packages(c("vegan", "psych", "lme4", "DHARMa", "performance", "MuMIn", "lmtest", "car", "robustlmm", "ggeffects"))

citation()
R.version.string

# load in each excel spreadsheet using file -> import dataset -> from excel


# PILOT SURVEYS species accumulation curves #
library(vegan)
# reshaping quadrat x species matrix
matrix_morning <- table(pilot_morningside$`Quadrat ID`, pilot_morningside$Species)
matrix_oldcal <- table(pilot_new_calton$`Quadrat ID`, pilot_new_calton$Species)
matrix_newcal <- table(pilot_old_calton$`Quadrat ID`, pilot_old_calton$Species)
# species accumulation curve 
for_plot_morn <- specaccum(matrix_morning, method = "random", permutations = 999)
for_plot_oldcal <- specaccum(matrix_oldcal, method = "random", permutations = 999)
for_plot_newcal <- specaccum(matrix_newcal, method = "random", permutations = 999)
plot(for_plot_morn, xlab = "Quadrats", ylab = "Species richness", main = "Morningside", cex.lab = 1.2)
plot(for_plot_oldcal, xlab = "Quadrats", ylab = "Species richness", main = "Old Calton", cex.lab = 1.2)
plot(for_plot_newcal, xlab = "Quadrats", ylab = "Species richness", main = "New Calton", cex.lab = 1.2)



# calculating mean and sd in vegetation height to assign mean vegetation height and sd to each burial ground, with 4 measurements in each quadrat averaged first and then quadrat averages averaged again for all quadrats within each site
vegheights <- vegetation_height_leaf_blade
vegheights$quadrat_mean <- rowMeans(vegheights[, c("Vegetation Height 1","Vegetation Height 2","Vegetation Height 3","Vegetation Height 4")], na.rm=TRUE)
site_mean <- tapply(vegheights$quadrat_mean, vegheights$Site, mean)
site_sd <- tapply(vegheights$quadrat_mean, vegheights$Site, sd)
print(site_mean)
print(site_sd)
# these were manually input into management_proxies datasheet
# proportion of quadrats with clipped leaf blades was also manually calculated

# adding column specifying survey round into species x percentage cover by quadrat by site for first and second survey round data
first_round_data$survey_round <- "first"
second_round_data$survey_round <- "second"
# combined dataframe
all_data <- rbind(first_round_data, second_round_data)

### SPECIES RICHNESS ###
# richness per site ##
richness_first <- tapply(first_round_data$Species, first_round_data$Site, function(x) length(unique(na.omit(x))))
richness_second <- tapply(second_round_data$Species, second_round_data$Site, function(x) length(unique(na.omit(x))))
# into dataframe
richness_df_sites <- data.frame(site = c(names(richness_first), names(richness_second)), richness = c(richness_first, richness_second), round = rep(c("First Round", "Second Round"), c(length(richness_first), length(richness_second))))

### SHANNON DIVERSITY ###
# making species x cemetery matrix 
species_matrix_first <- xtabs(`Percentage Cover` ~ Site + Species, data = first_round_data)
species_matrix_second <- xtabs(`Percentage Cover` ~ Site + Species, data = second_round_data)
# calculating shannon diversity
library(vegan)
# first round shannon diversity
shannon_first <- diversity(species_matrix_first, index = 'shannon')
# second round shannon diversity
shannon_second <- diversity(species_matrix_second, index = 'shannon')
# into data frame
shannon_df_sites <- data.frame(site = c(names(shannon_first), names(shannon_second)), shannon_diversity = c(shannon_first, shannon_second), round = rep(c("First Round", "Second Round"), c(length(shannon_first), length(shannon_second))))

### PROPORTION OF NON NATIVE SPECIES ###
# (lookup table created in excel with each unique species identified across all burial grounds during both survey rounds, with native status manually assigned referencing the PLANTATT database (Hill et al., 2004))
# assigning native status to each species record for both survey rounds
first_round_data <- merge(first_round_data, lookup[, c("Species", "NS")], by.x = "Species", by.y = "Species", all.x = TRUE)
second_round_data <- merge(second_round_data, lookup[, c("Species", "NS")], by.x = "Species", by.y = "Species", all.x = TRUE)
# NS.y column contains assigned status where N = native, AR = archaeophyte, AN = neophyte
# creating data frame for the species richness pool for each burial ground, with assigned native status
site_species_first <- unique(first_round_data[, c("Site", "Species", "NS")])
site_species_second <- unique(second_round_data[, c("Site", "Species", "NS")])
# calculating proportion of non-native species as number of unique AR and AN classified species compared to total number of unique species for each burial ground per survey round
proportion_nonative_first <- with(site_species_first, tapply(NS != "N", Site, mean, na.rm =TRUE))
proportion_nonative_second <- with(site_species_second, tapply(NS != "N", Site, mean, na.rm =TRUE))
# merging the prop. of non-native species from both survey rounds for all burial grounds by survey round
non_native_df <- data.frame(site = c(names(proportion_nonative_first), names(proportion_nonative_second)), proportion_non_native = c(proportion_nonative_first, proportion_nonative_second), round = rep(c("First Round", "Second Round"), c(length(proportion_nonative_first), length(proportion_nonative_second))))

### Creating data frame for GLMMs ###
# merging site x richness, site x shannon diversity, and site x proportion of non native species, by survey round.
glmm_data <- merge(richness_df_sites, shannon_df_sites, by = c("site", "round"))
glmm_data <- merge(glmm_data, non_native_df, by = c("site", "round"))
# merging site level data gathered, loaded in from excel with one site level observation per burial ground for both survey rounds
glmm_data <- merge(glmm_data, site_level_data, by.x = "site", by.y = "site_name")
# creating column for julian day for each of the two repeated observations for biodiversity metrics calculated
glmm_data$julian_day <- ifelse(glmm_data$round == "First Round", glmm_data$first_survey_julian, glmm_data$second_survey_julian)

### MANAGEMENT INTENSITY PCA ###
# removing site column
management_proxies <- management_proxies[, -1]
# kaiser meyer olkin test - sampling adequacy
library(psych)
KMO(management_proxies)
# bartletts test for sphericity - suitable correlation between variables
cortest.bartlett(cor(management_proxies), n = 30)
# PCA with each proxy scaled and centred
management_pca <- prcomp(management_proxies, center = TRUE, scale.= TRUE)
summary(management_pca)
# checking PC1 explains sufficient variance
screeplot(management_pca, type = "lines", main = "", xlab = "Principal Component", ylab = "Eigenvalue")
# checking the pc1 loading for each variable
round(management_pca$rotation[,1], 2)
# applying pca score to each burial ground then merging to glmm data 
site_level_data$management_PC1 <- management_pca$x[, 1]
glmm_data <- merge(glmm_data, site_level_data[, c("site_name", "management_PC1")], by.x = "site", by.y = "site_name", all.x = TRUE)


### SCALING ALL PREDICTORS ##
# calculating burial ground age from year of establishment
glmm_data$final_age <- 2026 - glmm_data$age
# z-scoring predictor values for each burial ground
glmm_data$size_z <- as.numeric(scale(glmm_data$size))
glmm_data$canopy_z <- as.numeric(scale(glmm_data$tree_canopy_percentage_cover))
glmm_data$julian_z <- as.numeric(scale(glmm_data$julian_day))
glmm_data$final_age_z <- as.numeric(scale(glmm_data$final_age))

### REFORMATTING PROPORTION NON-NATIVE SPECIES FOR GLMMS ###
# creating status column to assign species as 'native' (classified as N) or 'non-native' (classfied as AR and AN)
site_species_first$status <- ifelse(site_species_first$NS == "N", "native", "non_native")
site_species_second$status <- ifelse(site_species_second$NS == "N", "native", "non_native")
site_species_first$survey_round <- "First Round"
site_species_second$survey_round <- "Second Round"
# aggregating first round and second round unique species pools for each burial ground into a single data frame 
species_counts <- aggregate(Species ~ Site + survey_round + status, data = rbind(site_species_first, site_species_second), FUN = length)
# reformatting so unique non-native and native species counts are separate columns, by site and survey round
species_counts <- reshape(species_counts, idvar = c("Site", "survey_round"), timevar = "status", direction = "wide")
# NA for where non native species counts are zero, so replacing the NA with 0
species_counts$Species.non_native[is.na(species_counts$Species.non_native)] <- 0
# merging into dataframe for GLMMs 
glmm_data <- merge(glmm_data, species_counts, by.x = c("site", "round"), by.y = c("Site", "survey_round"), all.x = TRUE)
names(glmm_data)[names(glmm_data) == "Species.native"] <- "native_species"
names(glmm_data)[names(glmm_data) == "Species.non_native"] <- "non_native_species"
# creating binomial response object 
binomial_response <- cbind(glmm_data$non_native_species, glmm_data$native_species)


### CHECKING CORRELATIONS BETWEEN PREDICTORS ###
glmm_data$type <- as.factor(glmm_data$type)
predictors <- glmm_data[!duplicated(glmm_data$site), c("management_PC1", "age", "size", "tree_canopy_percentage_cover", "habitat_types", "julian_day")]
predictor_correlation_matrix <- round(cor(predictors, use = "complete.obs", method = "pearson"), 2)
plot(glmm_data$habitat_types, glmm_data$tree_canopy_percentage_cover)

### GLMMS ######################################
library(lme4)
library(DHARMa)
library(performance)
library(MuMIn)
library(lmtest)
library(car)

# adding in log10 area then scaling due to species-area relationships typically following a non-linear logathrimic relationship with greater increases in species richness for increased size when sizes are relatively small 
glmm_data$log_area <- log10(glmm_data$size)
glmm_data$log_area_z <- scale(glmm_data$log_area)
# checking distribution of burial grounds based on their age and size, coloring by type as there is clear separation based on churchyards (old and small) and cemeteries (young and large)
plot(glmm_data$size, glmm_data$final_age, col = as.factor(glmm_data$type), pch = 19, xlab = "Cemetery Size", ylab = "Cemetery Age")
legend("topright", legend = levels(as.factor(glmm_data$type)), col = 1:length(levels(as.factor(glmm_data$type))), pch = 19)
# type (cemetery and churchyard) included as an explanatory predictor as capturing differences in both age and size
# these included as alternative predictors in AICc model selection 

### Species richness ###########
# base model 
m_rich <- glmer(richness ~ management_PC1 + final_age_z + size_z + canopy_z + julian_z + (1|site), data = glmm_data, family = 'poisson')
# model with log transformed burial ground size
m_rich_log_size <- glmer(richness ~ management_PC1 + final_age_z + log_area_z + canopy_z + julian_z + (1|site), data = glmm_data, family = 'poisson')
# AICc model selection
model.sel(m_rich, m_rich_log_size, rank = "AICc")
# diagnostic plots
check_model(m_rich)
check_overdispersion(m_rich)
# model output
summary(m_rich)
# 95% confidence interval
confint(m_rich)
# variance inflation factors for each predictor 
vif(m_rich)

### Shannon diversity ############
# candidate model analysis (same format as species richness)
# uses maximum likelihood instead of restricted maximum likelihood due to differences between predictors
m_shan <- lmer(shannon_diversity ~ management_PC1 + final_age_z + size_z + canopy_z + julian_z + (1|site), data = glmm_data, REML = FALSE)
m_shan_log_size <- lmer(shannon_diversity ~ management_PC1 + final_age_z + log_area_z + canopy_z + julian_z + (1|site), data = glmm_data, REML = FALSE)

model.sel(m_shan, m_shan_log_size, rank = "AICc")
# using restricted maximum likelihood for final model selected
m_shan <- lmer(shannon_diversity ~ management_PC1 + final_age_z + size_z + canopy_z + julian_z + (1|site), data = glmm_data, REML = TRUE)
# model diagnostics plot 
check_model(m_shan)
# checking assumptions for homoscedasticity and normality of residuals 
check_heteroscedasticity(m_shan)
check_normality(m_shan)
# assumptions violated so refitting using robust lmm to downweight atypical residuals
library(robustlmm)
m_shannon_robust <- rlmer(shannon_diversity ~ management_PC1 + final_age_z + size_z + canopy_z + julian_z + (1 | site), data = glmm_data)
# model output
summary(m_shannon_robust)
check_model(m_shannon_robust)
# 95% confidence interval
confint(m_shannon_robust)
# variance inflation factors for each predictor
vif(m_shan)

### Proportion of Non-Native Species ##############
# candidate model analysis
m_nonnative <- glmer(binomial_response ~ management_PC1 + final_age_z + size_z + canopy_z + julian_z + (1|site), family = binomial, data = glmm_data)
m_nonnative_log_size <- glmer(binomial_response ~ management_PC1 + final_age_z + log_area_z + canopy_z + julian_z + (1|site), family = binomial, data = glmm_data)
model.sel(m_nonnative, m_nonnative_log_size, rank = "AICc")
# model output
summary(m_nonnative)
# model diagnostics
# simulating residuals 1000 times to conduct diagnostics for uniformity assumptions and dispersion
residuals_non_native <- simulateResiduals(m_nonnative, n = 1000)
plot(residuals_non_native)
testUniformity(residuals_non_native)
testDispersion(residuals_non_native)
# 95% confidence interval
confint(m_nonnative)
# variance of inflation factors between each predictor
vif(m_nonnative)
check_model(m_nonnative)
exp(0.39)
1 - exp(-0.20)

### POWER ANALYSIS ################################################
# using species richness model #
# using my actual predictor data
power_model <- glmm_data
power_model$site <- factor(power_model$site)
# using actual site sd as a random effect from model
site_sd <- as.data.frame(VarCorr(m_rich))$sdcor[1]
# selected effect sizes for main predictors (manually set each predictor effect ranging from 0.05 to 0.2)
management_effect <- 0.2
age_effect <- 0.2
size_effect <- 0.2
# using actual effect size for co variates
canopy_effect <- fixef(m_rich)["canopy_z"]
julian_effect <- fixef(m_rich)["julian_z"]
# using actual log expected species richness when predictors are held at mean found from final model as baseline richness
baseline_richness <- 3.312622
# setting number of simulation 
nsim <- 1000
# storage vectors for p-values
p_management <- numeric(nsim)
p_age <- numeric(nsim)
p_size <- numeric(nsim)
# simulation
for(i in 1:nsim){
  site_effect <- rnorm(nlevels(power_model$site), mean = 0, sd = site_sd)
  site_random <- site_effect[as.integer(power_model$site)]
  total_richness <- baseline_richness + management_effect * power_model$management_PC1 + age_effect * power_model$final_age_z + size_effect * power_model$size_z + canopy_effect * power_model$canopy_z + julian_effect * power_model$julian_z + site_random
  lambda <- exp(total_richness)
  power_model$richness_simulated <- rpois(nrow(power_model), lambda)
  model <- try(glmer(richness_simulated ~ management_PC1 + final_age_z + size_z + canopy_z + julian_z + (1|site), family = poisson, data = power_model), silent = TRUE)
  model_coef <- summary(model)$coefficients
  p_management[i] <- model_coef["management_PC1", "Pr(>|z|)"]
  p_age[i] <- model_coef["final_age_z", "Pr(>|z|)"]
  p_size[i] <- model_coef["size_z", "Pr(>|z|)"]
}

mean(p_management < 0.05, na.rm = TRUE)
mean(p_age < 0.05, na.rm = TRUE)
mean(p_size < 0.05, na.rm = TRUE)


### COMMUNITY COMPOSTION ANALYSIS ###################################

# creating sample id column to encode cemetery site and then survey round
all_data$site_survey_id <- paste(all_data$Site, all_data$survey_round, sep = "__")
# creating community matrix
community_matrix <- xtabs(`Percentage Cover` ~ site_survey_id + Species, data = all_data)
# applying square root transformation to dampen influence of very high percentage covers but retains quantitative differences in species abundance among samples.
community_matrix_sqrt <- sqrt(community_matrix)
# calculating bray-curtis dissimilarity matrices for first round and second round data by site (but in a combined matrix at this point)
library(vegan)
bray_curtis_dissimilarity <- vegdist(community_matrix_sqrt, method = "bray")
print(bray_curtis_dissimilarity)
# conducting two dimensional NMDS ordination with 200 starts to find best-fitting solution to visualise bray-curtis dissimilarities
nmds <- metaMDS(bray_curtis_dissimilarity, k = 2, trymax = 200, trace = FALSE)
nmds_scores <- as.data.frame(scores(nmds, display = "sites"))
# creating matching column for site_survey_id in glmm_data dataframe
glmm_data$survey_round <- ifelse(glmm_data$round == "First Round", "first", "second")
glmm_data$site_survey_id <- paste(glmm_data$site, glmm_data$survey_round, sep = "__")
# matching the scores to the site level data rows from glmm data by site_survey_id into a dataframe with NDMS scores
data_nmds <- glmm_data[match(rownames(nmds_scores), glmm_data$site_survey_id),]
identical(rownames(nmds_scores), data_nmds$site_survey_id)
#adding NMDS scores for both dimensions to each burial ground per survey round
data_nmds$NMDS1 <- nmds_scores$NMDS1
data_nmds$NMDS2 <- nmds_scores$NMDS2
#extracting into a first round dataframe and second round dataframe for PERMANOVA and NMDS ordination visualisation
data_first <- data_nmds[data_nmds$survey_round == "first", ]
data_second <- data_nmds[data_nmds$survey_round == "second", ]
#extracting bray curtis dissimilarities into first round dataframe and second round dataframe
bray_matrix <- as.matrix(bray_curtis_dissimilarity)
bray_first <- as.dist(bray_matrix[data_nmds$survey_round == "first", data_nmds$survey_round == "first"])
bray_second <- as.dist(bray_matrix[data_nmds$survey_round == "second", data_nmds$survey_round == "second"])
# conducting permanova separately for first round and second round, with 999 permutations with marginal calculations which tests associations with each predictor after accounting for all other predictors 
permanova_first <- adonis2(bray_first ~ management_PC1 + size_z + final_age_z + canopy_z, data = data_first, permutations = 999, by = "margin")
permanova_second <- adonis2(bray_second ~ management_PC1 + size_z + final_age_z + canopy_z, data = data_second, permutations = 999, by = "margin")
# conducting two dimensional NMDS ordination with 200 starts to find best-fitting solution to visualise bray-curtis dissimilarities for each matrix by survey round 
nmds_first <- metaMDS(bray_first, k = 2, trymax = 200, trace = FALSE)
nmds_second <- metaMDS(bray_second, k = 2, trymax = 200, trace = FALSE)
# calculating NMDS ordination stress and plotting stress to determine whether the two dimensional NMDS accurately captures 
nmds_first$stress
nmds_second$stress
stressplot(nmds_first)
stressplot(nmds_second)
# doing envfit() permutation test to visualize significant environmental variables to NMDS ordination plots
env_first <- envfit(nmds_first, data_first[, c("management_PC1", "size_z", "final_age_z", "canopy_z")], permutations = 999)
env_second <- envfit(nmds_second, data_second[, c("management_PC1", "size_z", "final_age_z", "canopy_z")], permutations = 999)
env_first
env_second
# visualising NMDS ordination with fitted significant environmental vectors
# function for applying gradient to burial ground site points based on their predictor values
plot_gradient <- function(nmds_object, variable, title, arrow_label) {
  cols <- hcl.colors(100, "Viridis")[cut(variable, 100, labels = FALSE)]
  plot(nmds_object, type = "n", main = title)
  points(nmds_object, display = "sites", pch = 19, col = cols, cex = 1.1)
  fitted <- envfit(nmds_object, setNames(data.frame(variable), arrow_label), permutations = 999)
  rownames(fitted$vectors$arrows) <- arrow_label
  plot(fitted, add = TRUE, p.max = 0.05)
}
# plotting NMDS ordinations for each predictor 
# first round
plot_gradient(nmds_first, data_first$management_PC1, "First round: Management Intensity", "")
legend("topright", legend = c("Low", "High"), col = hcl.colors(100, "Viridis")[c(1, 100)], pch = 19, bty = "n")
plot_gradient(nmds_first, data_first$size_z, "First round: Size", "Size")
legend("topright", legend = c("Small", "Large"), col = hcl.colors(100, "Viridis")[c(1, 100)], pch = 19, bty = "n")
plot_gradient(nmds_first, data_first$final_age_z, "First round: Age", "Age")
legend("topright", legend = c("Young", "Old"), col = hcl.colors(100, "Viridis")[c(1, 100)], pch = 19, bty = "n")
plot_gradient(nmds_first, data_first$canopy_z, "First round: Percentage canopy cover", "Canopy cover")
legend("topright", legend = c("Low", "High"), col = hcl.colors(100, "Viridis")[c(1, 100)], pch = 19, bty = "n")
# second round 
plot_gradient(nmds_second, data_second$management_PC1, "Second round: Management Intensity", "")
legend("topright", legend = c("Low", "High"), col = hcl.colors(100, "Viridis")[c(1, 100)], pch = 19, bty = "n")
plot_gradient(nmds_second, data_second$size_z, "Second round: Size", "Size")
legend("topright", legend = c("Small", "Large"), col = hcl.colors(100, "Viridis")[c(1, 100)], pch = 19, bty = "n")
plot_gradient(nmds_second, data_second$final_age_z, "Second round: Age", "Age")
legend("topright", legend = c("Young", "Old"), col = hcl.colors(100, "Viridis")[c(1, 100)], pch = 19, bty = "n")
plot_gradient(nmds_second, data_second$canopy_z, "Second round: Percentage Canopy Cover", "Canopy cover")
legend("topright", legend = c("Low", "High"), col = hcl.colors(100, "Viridis")[c(1, 100)], pch = 19, bty = "n")


#### RANK-ABUNDANCE CURVE ########################################################################################

# summing percentage cover found for each species across all quadrats for all burial grounda across both survey rounds
taxon_cover <- tapply(all_data$`Percentage Cover`, all_data$Species, sum, na.rm = TRUE)
# removing missing/zero totals because there was one quadrat in two burial grounds that were empty so had no percentage cover values
taxon_cover <- taxon_cover[!is.na(taxon_cover) & taxon_cover > 0]
# ranking species from hight to lowest total percentage cover
taxon_cover <- sort(taxon_cover, decreasing = TRUE)
# calculating the percentage contribution of each species total percentage cover across surveys to the total recorded percentage cover across all species 
relative_cover <- taxon_cover / sum(taxon_cover) * 100
# ranking species from highest relative cover contribution to lowest
taxon_rank <- seq_along(relative_cover)

# plotting rank-abundance curve with x axis as the ranked relative contribution for each species and y axis as the relative % contribution to total percentage cover across all species
plot(taxon_rank, relative_cover, type = "o", pch = 16, cex = 0.45, lwd = 1, xlab = "Taxon rank", ylab = "Contribution to total recorded cover (%)", main = "")
# labelling with species name for the 10 most abundance species on the rank-abundance curve
text(taxon_rank[1:10], relative_cover[1:10], labels = names(relative_cover)[1:10], pos = c(4, 4, 4, 4, 4, 4, 4, 4, 4, 4), offset = 0.5, cex = 0.6, font = 3)
# total percent contribution to the total percentage cover across all species calculated by summing the 10 dominant species
sum(relative_cover[1:10])

# complimentary finding the frequency in which the 10 most dominant taxa in terms of percentage cover contribution occur across burial grounds
top_10_taxa <- names(relative_cover)[1:10]
table(unique(all_data[all_data$Species %in% top_10_taxa, c("Site", "Species")])$Species)[top_10_taxa]



###### DESCRIPTIVE STATS ##########################################################################################
summary(glmm_data$final_age)
summary(glmm_data$size)
summary(glmm_data$management_PC1)
summary(glmm_data$tree_canopy_percentage_cover)
sd(glmm_data$tree_canopy_percentage_cover)

summary(glmm_data$richness)
sd(glmm_data$richness)
summary(glmm_data$native_species)
summary(glmm_data$non_native_species)

sum(is.na(lookup$NS))
sum(lookup$NS == "N", na.rm = TRUE)
sum(lookup$NS %in% c("AR", "AN"), na.rm = TRUE)

summary(glmm_data$shannon_diversity)
sd(glmm_data$shannon_diversity, na.rm = TRUE)



### finding sum of burial ground size ###
sum(glmm_data$size)
# 1594276


###### BETA DISPERSION TEST #########################################################################################

# as management intensity needs to be a discrete grouping factor, assigning low management and high management based on whether each burial ground score is > PC1 median or < PC1 median :)
site_level_data$management_group <- ifelse(site_level_data$management_PC1 >= median(site_level_data$management_PC1), "High", "Low")
# checking site names match and in same order
labels(bray_first)
site_level_data$site_name

# beta dispersion test for first round bray curtis dissimilarity matrix, then permutation test with the distances between averaged to determine whether bray-curtis dissimilarities converge under high management intensity
bd_first <- betadisper(bray_first, site_level_data$management_group)
permutest(bd_first, permutations = 999)
aggregate(bd_first$distances, by = list(site_level_data$management_group), FUN = mean)
# same for second round 
bd_second <- betadisper(bray_second, site_level_data$management_group)
permutest(bd_second, permutations = 999)
aggregate(bd_second$distances, by = list(site_level_data$management_group), FUN = mean)