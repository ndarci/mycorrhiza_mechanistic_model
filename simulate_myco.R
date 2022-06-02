# implement EEB 219B mycorrhiza model

library(tidyverse)
library(deSolve)
library(grid)
library(gridExtra)
library(cowplot)
theme_set(theme_bw())

setwd('~/src/eeb219b_mycorrhiza_model/')

# define function to update derivatives at each time step
myco_update <- function(tvec, sv, params) {
  # update SVs from last time step
  for (i in 1:length(sv)) {
    assign(names(sv)[i], sv[i])
  }
  # remember parameter values
  for (i in 1:length(params)) {
    assign(names(params)[i], params[i])
  }
  
  # define time derivatives for each variable
  dCT.dt = a - k*CT - M*CT*NTmin/NT
  dNF.dt = g - k*NF - M*NF*CFmin/CF
  dNT.dt = m - k*NT + M*NF*CFmin/CF
  dCF.dt = b - k*CF + M*CT*NTmin/NT
  dM.dt = M*(1-M)*CFmin/CF - d*M
  
  return(list(c(dCT = dCT.dt,
              dNT = dNT.dt,
              dCF = dCF.dt,
              dNF = dNF.dt,
              dM = dM.dt)))
}

# function to simulate the system and return dataframe of results
run_sim <- function(tvec, svs, params) {
  # simulate the model over time
  res = lsoda(sv, tvec, myco_update, params) %>% data.frame(.)
  
  # augment with important flow amounts
  res['CT_to_tree'] = params['k'] * res$CT
  res['CT_to_myco'] = res$M * res$CT * params['NTmin'] / res$NT
  res['NF_to_fungus'] = params['k'] * res$NF
  res['NF_to_myco'] = res$M * res$NF * params['CFmin'] / res$CF
  
  # rearrange to get nutrient dynamics over time
  nut_time = res %>% pivot_longer(cols = c('CT', 'NT', 'CF', 'NF',
                                           'CT_to_tree', 'CT_to_myco',
                                           'NF_to_fungus', 'NF_to_myco'),
                             names_to = 'variable',
                             values_to = 'g_per_L') %>% 
  mutate(nutrient = substr(variable, 1, 1),
           species = substr(variable, 2, 2))
  
  return(list(res = res, nut_time = nut_time))
}

# define function to get steady state values for summary-level stats of one sim
get_steady_state_vals <- function(rawres) {
  NSTEPS = 2
  # get last few time steps, ignoring time column
  lastNsteps = tail(rawres[,2:ncol(rawres)], NSTEPS)
  # find the final value of each SV/flow, if it was only slightly different from the previous time step
  diff = lastNsteps[nrow(lastNsteps),] - lastNsteps[1,]
  ss = as.numeric(lastNsteps[nrow(lastNsteps),])
  ss[diff > 0.001] <- NA
  names(ss) = names(lastNsteps)
  return(ss)
}

# define function to plot time series of nutrient state variables
plot_SVs_time <- function(timedata, params) {
  timedata = timedata[timedata$variable %in% c('CT', 'CF', 'NT', 'NF'),]
  title = ''
  if(params['a'] > 1) {
    title = paste0(title, 'Abundant C, ')
  } else {
    title = paste0(title, 'Scarce C, ')
  }
  if(params['g'] > 1) {
    title = paste0(title, 'Abundant N')
  } else {
    title = paste0(title, 'Scarce N')
  }
  nutplot = ggplot(timedata, aes(x = time, y = g_per_L, 
                                 color = species)) + 
    scale_color_manual(breaks = c('T', 'F'),
                       values = c('forestgreen', 'lightcoral'),
                       labels = c('Tree', 'Fungus'),
                       name = 'Species') +
    geom_line(size = 1) +
    facet_grid(nutrient~., scales = 'free_y', 
               labeller = labeller(nutrient = c(C = 'Free C', N = 'Free N'))) +
    ggtitle(title) + ylab('g/L') +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlim(0, 12)
  return(nutplot)
}

# plot flows from crucial compartments to self vs. myco
plot_flows_time <- function(timedata, params) {
  timedata = timedata[timedata$variable %in% c('CT_to_tree', 'CT_to_myco',
                                               'NF_to_fungus', 'NF_to_myco'),]
  timedata['variable'] = recode(timedata$variable, CT_to_tree = 'to_self',
                                                    NF_to_fungus = 'to_self',
                                                    CT_to_myco = 'to_myco',
                                                    NF_to_myco = 'to_myco')
  timedata['colorcode'] = timedata$species
  timedata[timedata$variable == 'to_myco', 'colorcode'] <- 'M'
  flowplot = ggplot(timedata, aes(x = time, y = g_per_L,
                                 color = colorcode)) +
    scale_color_manual(name = '',
                       breaks = c('T', 'F', 'M'),
                       values = c('forestgreen', 'lightcoral', 'dodgerblue1'),
                       labels = c('Tree C to self', 'Fungus N to self', 'To mycorrhiza')) +
    geom_line(size = 1) +
    facet_grid(species~., scales = 'free_y') +
    ggtitle(paste0('m=', params['m'],
                   ', a=', params['a'],
                   ', b=', params['b'],
                   ', g=', params['g'])) +
    xlim(0, 12)
  return(flowplot)
}

# function to plot time series of mycorrhiza state variable
plot_myco_time <- function(rawres) {
  # plot M over time
  mplot = ggplot(rawres, aes(x = time, y = M)) + 
    geom_line(color = 'dodgerblue1', size = 1) +
    xlim(0, 12) + 
    xlab('Time') + ylab('% of roots colonized')
  return(mplot)
}

# function to plot dynamics of M with changes in nutrient SVs
plot_myco_vs_nutrients <- function(rawres, params) {
  title = ''
  if(params['a'] > 1) {
    title = paste0(title, 'Abundant C, ')
  } else {
    title = paste0(title, 'Scarce C, ')
  }
  if(params['g'] > 1) {
    title = paste0(title, 'Abundant N')
  } else {
    title = paste0(title, 'Scarce N')
  }
  # plot CT vs NT vs M
  tree_vs_m = ggplot(rawres, aes(x = CT, y = NT, color = M)) + geom_point(size = 3) +
    xlab('Tree free C') + ylab('Tree free N') +
    scale_color_continuous(name = '% roots colonized', limits = c(0, 1))
  fungus_vs_m = ggplot(rawres, aes(x = CF, y = NF, color = M)) + geom_point(size = 3) +
    xlab('Fungus free C') + ylab('Fungus free N') + 
    scale_color_continuous(name = '% roots colonized', limits = c(0, 1))
  mleg = cowplot::get_legend(fungus_vs_m)
  return(list(plot = grid.arrange(tree_vs_m + theme(legend.position = 'none'), 
                                  fungus_vs_m + theme(legend.position = 'none'), 
                                  ncol = 2,
                                  top = textGrob(title)), 
              legend = mleg))
}

# define time steps for the simulation
tvec <- seq(0, 100, by = 0.1)

# define state variables (SVs) and their initial values
CT0 = 1
NT0 = 1
CF0 = 1
NF0 = 1
M0 = 0.1
sv = c(CT = CT0, NT = NT0, 
       CF = CF0, NF = NF0,
       M = M0)

# iterate over a series of the most important inputs
avec = seq(0, 100, 5)
gvec = seq(0, 100, 5)
avec[1] = 1
gvec[1] = 1
# keep track of inputs/outputs
ionames = c('alpha', 'gamma', paste0('ss.', c('CT', 'NT', 'CF', 'NF', 'M',
            'CT_to_tree', 'CT_to_myco', 'NF_to_fungus', 'NF_to_myco')))
io = data.frame(matrix(nrow = length(avec)*length(gvec), ncol = length(ionames)))
names(io) = ionames
ioiter = 1
# keep track of the most extreme case plots
alltimeseries = vector('list', 4)
allphase = vector('list', 4)
plotiter = 1
for (thisa in avec) {
  for (thisg in gvec) {
    # define parameters and their values on this iteration
    params = c(NTmin = 1, # minimum free [N] required for healthy tree
               CFmin = 1, # minimum free [C] required for healthy fungus
               a = thisa, # C input to tree per time
               b = 0, # C input to fungus per time
               g = thisg, # N input to fungus per time
               m = 0, # N input to tree per time
               k = 1, # rate that tree and fungus take nutrients for themselves
               d = 0.5) # rate that myco connections wear off over time
    
    # run the simulation with these parameters
    sim_result = run_sim(tvec, svs, params)
    res = sim_result$res
    nut_time = sim_result$nut_time
    
    # get steady state values
    steady_sv = get_steady_state_vals(res)
    
    # add to in/out dataframe
    io[ioiter,] = c(thisa, thisg, steady_sv)
    ioiter = ioiter+1
    
    # in the most extreme cases...
    if ((thisa == avec[1] | thisa == avec[length(avec)]) & 
      (thisg == gvec[1] | thisg == gvec[length(gvec)])) {
        # plot nutrients and M over time, then put the plots together
        nutplot = plot_SVs_time(nut_time, params)
        mplot = plot_myco_time(res)
        mycoplot = grid.arrange(nutplot + theme(legend.position = 'none',
                                          axis.title.x = element_blank()), 
                                mplot, 
                                nrow = 2, heights = c(2, 1))
        # save this plot to the plot vector
        alltimeseries[[plotiter]] = mycoplot
        
        # add a legend at the end of the plot vector
        if(thisa == avec[length(avec)] & thisg == gvec[length(gvec)]) {
          alltimeseries = append(list(cowplot::get_legend(nutplot)), alltimeseries)
        }
        
        # plot phase dynamics of nutrients
        phaseplotobj = plot_myco_vs_nutrients(res, params)
        allphase[[plotiter]] = phaseplotobj$plot
        
        # add legend at the end
        if(thisa == avec[length(avec)] & thisg == gvec[length(gvec)]) {
          allphase = append(list(phaseplotobj$legend), allphase)
        }
        
        plotiter = plotiter+1
      }
  }
}

# define mutualism/parasitism/commensalism based on outputs
io['tree_healthy'] = io$ss.CT > sv['CT'] & io$ss.NT > sv['NT']
io['fungus_healthy'] = io$ss.CF > sv['CF'] & io$ss.NF > sv['NF']
io = io %>% mutate(relationship = ifelse(tree_healthy & fungus_healthy, '+/+',
                                         ifelse(tree_healthy & !fungus_healthy, '+/-',
                                                ifelse(!tree_healthy & fungus_healthy, '-/+',
                                                'Both fail'))))

# plot steady state outputs over different inputs
relationshipplot = ggplot(io, aes(x = factor(alpha), y = factor(gamma), fill = relationship)) +
  geom_tile() +
  coord_fixed() + 
  scale_fill_manual(name = 'Relationship (T/F)',
                      breaks = c('+/+', '+/-', '-/+', 'Both fail'),
                      values = c('goldenrod1', 'mediumseagreen', 'skyblue3', 'thistle1')) +
  xlab("C input to tree per time (g/L)") + ylab("N input to fungus per time (g/L)") +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        legend.position = 'right',
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())
ggsave('relationship_heatmap.png', relationshipplot, dpi = 300)

m_heatmap = ggplot(io, aes(x = factor(alpha), y = factor(gamma), fill = ss.M)) +
  geom_tile() +
  coord_fixed() + 
  xlab("C input to tree per time (g/L)") + ylab("N input to fungus per time (g/L)") +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        legend.position = 'right',
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())
  
# put all the extreme-case nutrient plots together
coollayout = matrix(2:(length(alltimeseries)),
                    ncol = floor(sqrt(length(alltimeseries)-1)),
                    byrow = T)
coollayout = cbind(coollayout, rep(1, dim(coollayout)[1]))
coolwidths = c(rep(5, (dim(coollayout)[2] - 1)), 1)
multi = grid.arrange(grobs = alltimeseries,
             layout_matrix = coollayout,
             widths = coolwidths)
ggsave('extremes_myco_dynamics.png', multi, height = 9, width = 12, units = 'in', dpi = 300)

multi_phase = grid.arrange(grobs = allphase,
                           layout_matrix = coollayout,
                           widths = coolwidths)
ggsave('extremes_M_vs_nutrients.png', multi_phase, height = 9, width = 15, units = 'in', dpi = 300)





