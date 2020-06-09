# metadata
# - code associated with article: 
# -- "Black, Asian and Minority Ethnic groups in England are at increased risk of death from COVID-19: indirect standardisation of NHS mortality data"
# -- https://doi.org/10.12688/wellcomeopenres.15922.1
# - this code is written for R version 3.6.3
# - author - Dan Lewer (d.lewer@ucl.ac.uk)
# - version 2 (after peer review): 9 June 2020
# - sources:
# -- Covid-19 deaths data: https://www.england.nhs.uk/statistics/statistical-work-areas/covid-19-daily-deaths/
# -- census 2011 population estimates: table DC2101EW - Ethnic group by sex by age. nomisweb.co.uk

library(data.table)
library(RColorBrewer)
library(colorspace) # for 'lighten' function

#----------
# read data
#---------- 

pop <- fread("https://raw.githubusercontent.com/danlewer/covid19_ethnicity/master/nomis_reformatted.csv")
d.eth <- fread("https://raw.githubusercontent.com/danlewer/covid19_ethnicity/master/ethnicity_26apr2020.csv")
d.age <- fread("https://raw.githubusercontent.com/danlewer/covid19_ethnicity/master/age_26apr2020.csv")
d.region <- fread("https://raw.githubusercontent.com/danlewer/covid19_ethnicity/master/region_26apr2020.csv")
lookup_eth <- fread("https://raw.githubusercontent.com/danlewer/covid19_ethnicity/master/lookup_ethnicity.csv")
lookup_age <- fread("https://raw.githubusercontent.com/danlewer/covid19_ethnicity/master/lookup_age.csv")
lookup_region <- fread("https://raw.githubusercontent.com/danlewer/covid19_ethnicity/master/lookup_region.csv")

#---------------------------------------
# align number of deaths in each dataset
#---------------------------------------

eth_factor <- sum(d.eth$d) / sum(d.eth$d[!(d.eth$ethnicity %in% c('Not stated', 'No match'))])
d.eth[, d := d * eth_factor]
d.eth <- d.eth[!(ethnicity %in% c('Not stated', 'No match'))]

age_factor <- sum(d.age$d) / sum(d.eth$d)
d.age[, d := d / age_factor]

region_factor <- sum(d.region$d) / sum(d.eth$d)
d.region[, d := d / region_factor]

#---------------------------------------------------
# align categories in population data with NHSE data
#---------------------------------------------------

pop <- melt(pop, id.vars = c('ethnicity', 'region'), variable.name = 'age', value.name = 'p', variable.factor = F)
pop <- lookup_age[pop, on = 'age']
pop <- lookup_region[pop, on = 'region']
pop <- lookup_eth[pop, on = 'ethnicity']
pop <- pop[!is.na(age_nhse) & !is.na(region_nhse) & !is.na(ethnicity_nhse), .(p = sum(p)), c('age_nhse', 'region_nhse', 'ethnicity_nhse')]
setnames(pop, c('age_nhse', 'region_nhse', 'ethnicity_nhse'), c('age', 'region', 'ethnicity'))

#-----------------------------------------------------
# expected deaths - adjusted for age only (not region)
#-----------------------------------------------------

smr1 <- pop[, .(p = sum(p)), age][d.age, on = 'age']
smr1[, mr := d / p]
smr1 <- pop[, .(p = sum(p)), c('age', 'ethnicity')][smr1[, c('age', 'mr')], on = 'age']
smr1[, exp1 := mr * p]
smr1 <- smr1[, .(exp1 = sum(exp1)), ethnicity]

#--------------------------------
# estimate deaths by region x age
#--------------------------------

# assume that deaths have the same age breakdown within each region

d.age[, pc := d / sum(d)]
d.age[, d := NULL]
d.er <- CJ(age = d.age$age, region = d.region$region)
d.er <- d.age[d.er, on = 'age']
d.er <- d.region[d.er, on = 'region']
d.er[, est_d := d * pc]

#-------------------------------
# mortality rate by region x age
#-------------------------------

pop_ra <- pop[, .(p = sum(p)), c('region', 'age')]
pop_ra <- pop_ra[d.er, on = c('region', 'age')]
pop_ra[, mr := est_d / p] # higher rates in London by age

#----------------------------------------
# expected deaths by ethnicity and region
#----------------------------------------

pop <- pop_ra[, c('region', 'age', 'mr')][pop, on = c('age', 'region')]
pop[, exp := mr * p]
smr2 <- pop[, .(exp2 = sum(exp)), ethnicity]

#------------------------------------
# compare SMRs and summarise by group
#------------------------------------

smr <- d.eth[smr1, on = 'ethnicity'][smr2, on = 'ethnicity']
smr[, d := round(d, 0)]

short_names <- data.table(ethnicity = lookup_eth$ethnicity_nhse, short_name = lookup_eth$short_name)
short_names <- unique(short_names[!is.na(ethnicity)])

smr <- short_names[smr, on = 'ethnicity']

lookup_eth[, ethnicity := NULL]
setnames(lookup_eth, 'ethnicity_nhse', 'ethnicity')
lookup_eth <- unique(lookup_eth[complete.cases(lookup_eth)])
smr <- lookup_eth[, c('ethnicity', 'group')][smr, on = 'ethnicity']
smr <- rbind(smr, smr[, lapply(.SD, sum), group, .SDcols = c('d', 'exp1', 'exp2')], fill = T) # add totals
smr <- smr[order(group, d)]
smr[, ethnicity := ifelse(is.na(ethnicity), paste0('Total: ', group), ethnicity)]
bold <- is.na(smr$short_name)
smr[, short_names := ifelse(is.na(short_name), ethnicity, short_name)]
bold <- bold[smr$short_name != 'Other']
smr <- smr[short_names != 'Other']

vpt <- function(vx, vt) { # poisson test across vector
  pt <- function(x, t) {
    r <- poisson.test(x, t)
    c(x = x, t = t, rate = x/t, lower = r$conf.int[1], upper = r$conf.int[2], p_value = r$p.value)
  }
  t(mapply(pt, x = vx, t = vt))
}

ci1 <- vpt(smr$d, smr$exp1)[, c('rate', 'lower', 'upper', 'p_value')]
ci2 <- vpt(smr$d, smr$exp2)[, c('rate', 'lower', 'upper', 'p_value')]

#-----
# plot
#-----

cols1 <- brewer.pal(5, 'Set2')
cols2 <- lighten(cols1, amount = 0.8)
eg <- data.table(group = unique(smr$group), col1 = cols1, col2 = cols2)
smr <- eg[smr, on = 'group']

x1l <- 0:19
x1r <- 0:19 + 0.3
cix1 <- rowMeans(cbind(x1l, x1r))
x2l <- 0:19 + 0.3
x2r <- 0:19 + 0.6
cix2 <- rowMeans(cbind(x2l, x2r))
cil <- rowMeans(cbind(cix1, cix2))

ylims = range(rbind(ci1, ci2)[,-4])
yx <- c(0.5, 1, 2, 5, 10)

#png('covid_smr_ethnicity_27apr2020.png', height = 7.5, width = 10, units = 'in', res = 300)

par(xpd = NA, mar = c(10, 4, 2, 10))
plot(1, type = 'n', xlim = c(0, 22), ylim = log(ylims), axes = F, xlab = NA, ylab = NA)

rect(x1l, 0, x1r, log(ci1[,1]), col = smr$col1)
rect(x2l, 0, x2r, log(ci2[,1]), col = smr$col2)

arrows(cix1, log(ci1[,2]), cix1, log(ci1[,3]), length = 0.05, angle = 90, code = 3)
arrows(cix2, log(ci2[,2]), cix2, log(ci2[,3]), length = 0.05, angle = 90, code = 3)

axis(2, log(yx), yx, las = 2, pos = -1)
text(cil, log(ylims)[1], smr$short_names, srt = 90, adj = 1, font = ifelse(bold, 2, 1))
title(ylab = 'SMR')

segments(-1, 0, 21)

ly <- c(0, log(ylims[2]) * 0.5) 
ly <- seq(ly[1], ly[2], length.out = nrow(eg) + 1)

rect(22, ly[-length(ly)], 23, ly[-1], col = rev(eg$col1))
rect(23, ly[-length(ly)], 24, ly[-1], col = rev(eg$col2))
text(24.25, ly[-length(ly)] + diff(ly) / 2, rev(eg$group), adj = 0)
text(c(22.5, 23.5), max(ly) + 0.05, c('Adjusted for age', 'Adjusted for age and region'), srt = 90, adj = 0)

#dev.off()

#-----------
# make table
#-----------

fci <- function(x) {
  x <- format(round(x, 2), digits = 2, nsmall = 2)
  x <- paste0(x[,1], '(', x[,2], '-', x[,3], ')')
  x <- gsub(' ', '', x)
  gsub('\\(', ' \\(', x)
}
p_summarise <- function(p) {
  x <- round(p, 3)
  x[p < 0.001] <- '<0.001'
  x
}

nice_table <- smr[, c('group', 'short_names', 'd', 'exp1', 'exp2')]
nice_table[, smr1 := fci(ci1)]
nice_table[, smr2 := fci(ci2)]
nice_table[, d := format(d, big.mark = ',')]
nice_table[, exp1 := format(round(exp1, 1), digits = 1, nsmall = 1, big.mark = ',')]
nice_table[, exp2 := format(round(exp2, 1), digits = 1, nsmall = 1, big.mark = ',')]
nice_table[, x := 1]
nice_table[, x := cumsum(x), group]
nice_table[, x := x != 1]
nice_table$group[nice_table$x] <- ''
nice_table[, p1 := p_summarise(ci1[,4])]
nice_table[, p2 := p_summarise(ci2[,4])]
nice_table <- nice_table[, c('group', 'short_names', 'd', 'exp1', 'smr1', 'p1', 'exp2', 'smr2', 'p2')]
setnames(nice_table, c('group', 'short_names', 'd', 'exp1', 'smr1', 'p1', 'exp2', 'smr2', 'p2'),
         c('Group', 'Subgroup', 'Observed deaths', 'Expected deaths (adjusting for age)', 'SMR (adjusting for age)', 'p-value (adjusting for age)', 'Expected deaths (adjusting for age and region)', 'SMR (adjusting for age and region)', 'p-value (adjusting for age and region)'))

#fwrite(nice_table, 'ethnicity_smr_table_28apr2020.csv')
