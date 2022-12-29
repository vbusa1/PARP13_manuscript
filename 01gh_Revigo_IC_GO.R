# A plotting R script produced by the Revigo server at http://revigo.irb.hr/
# If you found Revigo useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

library( ggplot2 );
library(ggrepel)
library( scales );

# --------------------------------------------------------------------------

######### WT vs KO IC (IC2) ########
revigo.names <- c("term_ID","description","frequency","plot_X","plot_Y","log_size","value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0034976","response to endoplasmic reticulum stress",1.2175187418973,-5.72089191607853,-3.11136073326826,2.33645973384853,-9.99139982823808,0.814599010760252,0),
c("GO:1905897","regulation of response to endoplasmic reticulum stress",0.456569528211488,2.11202954329701,4.43514874229215,1.91381385238372,-6.09097914578884,0.826398848170483,0),
c("GO:0022900","electron transport chain",0.930049038949326,4.61460510289546,-3.29727937987509,2.22010808804005,-7.03198428600636,0.960689762850783,0.00512444),
c("GO:1900118","negative regulation of execution phase of apoptosis",0.124006538526577,-2.24331494619971,5.64907144109105,1.36172783601759,-5.02965312376991,0.748856719400529,0.02866141),
c("GO:0010941","regulation of cell death",8.87210416549236,5.32159203325518,1.48662445620919,3.19728055812562,-5.68613277963085,0.947182476711719,0.04432994),
c("GO:0006091","generation of precursor metabolites and energy",2.19829772842568,3.04533614577627,-6.13757502315202,2.59217675739587,-5.20411998265593,0.957928559133036,0.12986242),
c("GO:0006950","response to stress",18.5727974747759,-5.89570723050972,1.46667273121342,3.51798720302508,-5.08671609823958,0.920898880571446,0.14469001),
c("GO:0044281","small molecule metabolic process",8.73118764443943,0.508716436124404,-7.71372169943732,3.19033169817029,-5.80966830182971,0.961198943957789,0.14596989),
c("GO:0033554","cellular response to stress",8.38171467222817,-5.73244345089213,-1.76738890845029,3.17260293120986,-5.3269790928711,0.833810202207365,0.39805726),
c("GO:0030968","endoplasmic reticulum unfolded protein response",0.298743024632208,-5.55073722217786,-4.23878656458475,1.73239375982297,-9.53910215724345,0.743988957879842,0.44441253));

one.data <- data.frame(revigo.data);
names(one.data) <- revigo.names;
one.data <- one.data [(one.data$plot_X != "null" & one.data$plot_Y != "null"), ];
one.data$plot_X <- as.numeric( as.character(one.data$plot_X) );
one.data$plot_Y <- as.numeric( as.character(one.data$plot_Y) );
one.data$log_size <- as.numeric( as.character(one.data$log_size) );
one.data$value <- as.numeric( as.character(one.data$value) );
one.data$frequency <- as.numeric( as.character(one.data$frequency) );
one.data$uniqueness <- as.numeric( as.character(one.data$uniqueness) );
one.data$dispensability <- as.numeric( as.character(one.data$dispensability) );

ggplot(one.data, aes(x = -value, y = reorder(description, -value), 
                     color = log_size)) +
    geom_point(size = 5) +
    theme_classic() +
    scale_color_continuous("GO term\n frequency", limits = c(0,3.6)) +
    labs(y=NULL,x="-log(p)") +
    scale_x_continuous(lim = c(0, 10), breaks = seq(0, 10, 5)) +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          legend.text = element_text(size = 12))
ggsave("Revigo_IC2_WT.jpg", 
       width = 7, height = 3, 
       units = "in")



######### IMMUNE IC (IC3) ###########
revigo.names <- c("term_ID","description","frequency","plot_X","plot_Y","log_size","value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0002252","immune effector process",2.4012175187419,3.20843485999068,4.9024418058183,2.63042787502502,-42.4841261562883,0.843191005083917,0),
                     c("GO:0002376","immune system process",12.947409954343,-1.42277533262016,-8.41536881458904,3.36135002435227,-47.1985962899826,1,0),
                     c("GO:0009615","response to virus",1.88264472126712,-3.81590383289042,-6.73325930091213,2.52504480703685,-67.8894102897007,0.844182604779755,0),
                     c("GO:0016032","viral process",1.42043853221352,-6.92429418813962,2.98082349104779,2.40312052117582,-9.44249279809434,1,0),
                     c("GO:0032020","ISG15-protein conjugation",0.0338199650527028,4.60110371788809,-0.043392478333703,0.845098040014257,-7.54363396687096,0.991033036360388,0),
                     c("GO:0044419","biological process of interspecies interaction",8.57336114086015,0.852281910692247,-8.37187325467899,3.18241465243455,-10.7235381958268,1,0),
                     c("GO:0048525","negative regulation of viral process",0.507299475790542,5.77002984842702,2.98592410266736,1.95904139232109,-36.5114492834996,0.863657610801178,0),
                     c("GO:0050896","response to stimulus",44.9241869116735,1.43763219909277,-0.941129861511749,3.90151280912994,-28.5968794788242,1,0),
                     c("GO:0043903","regulation of biological process in symbiotic interaction",0.321289668000676,-0.222059903472274,-6.71988307988924,1.76342799356294,-30.9546770212133,0.966929538647436,0.02644702),
                     c("GO:0050792","regulation of viral process",0.884955752212389,-0.2253132104918,2.95550961037339,2.19865708695442,-31.7121982700698,0.963277133943988,0.02912362),
                     c("GO:0031347","regulation of defense response",3.4214531311651,3.90129185767291,-4.67730676475482,2.78390357927274,-15.9244530386075,0.780637289486942,0.03598432),
                     c("GO:0002697","regulation of immune effector process",1.96719463389888,6.98433051187871,-0.03134824899057,2.54406804435028,-11.7189666327523,0.787559520352563,0.03992004),
                     c("GO:0001817","regulation of cytokine production",4.0471224846401,-2.91927158274762,5.74451267537436,2.85672889038288,-13.8728952016352,0.893012360407774,0.04429517),
                     c("GO:0002682","regulation of immune system process",7.97587509159574,7.00194504020175,-1.84585335047645,3.15106325335375,-17.8927900303521,0.951535287508988,0.05083109),
                     c("GO:0048519","negative regulation of biological process",29.1471732145877,0.615614548964806,6.4716239797552,3.71365851620836,-9.43415218132648,0.940138567295417,0.07751109),
                     c("GO:0007259","receptor signaling pathway via JAK-STAT",0.253649737895271,-2.84254955117105,-0.227390596071504,1.66275783168157,-9.80687540164554,0.890022226771769,0.08859036),
                     c("GO:0048583","regulation of response to stimulus",21.7687841722564,2.2766024040861,1.74236595769284,3.58692470814482,-9.13786862068696,0.943134562108693,0.10573422),
                     c("GO:0034097","response to cytokine",4.43605208274618,-6.83497915778261,-1.8984638513851,2.89652621748956,-16.0264103765727,0.864344367523486,0.12435655),
                     c("GO:0045343","MHC class I biosynthetic process regulation",0.0394566258948199,-1.36021618817191,6.076150547109,0.903089986991944,-5.12436006299583,0.954520080200032,0.12557846),
                     c("GO:0060700","regulation of ribonuclease activity",0.045093286736937,-4.64526103978476,4.20756262389004,0.954242509439325,-8.44369749923271,0.947457805038007,0.12711527),
                     c("GO:1900118","negative regulation of apoptosis execution",0.124006538526577,6.29354838867469,1.97445327993324,1.36172783601759,-5.44490555142168,0.939132095587588,0.14586708),
                     c("GO:0006952","defense response",7.8631418747534,-5.50864789283496,-5.4430970369917,3.14488541828714,-65.3098039199715,0.885198000936471,0.1557734),
                     c("GO:0009607","response to biotic stimulus",7.98151175243786,-3.95203794946829,-5.05193374980629,3.15136985024746,-57.0301183562535,0.918276595905058,0.1738162),
                     c("GO:0009605","response to external stimulus",13.3532495349755,-6.07528487498141,-4.26275481399471,3.3747483460101,-41.9430951486635,0.911013506309087,0.19408471),
                     c("GO:0006950","response to stress",18.5727974747759,-5.33838675655729,-4.04989950727731,3.51798720302508,-31.9393021596464,0.905643307580732,0.23842427),
                     c("GO:0007166","cell surface receptor signaling pathway",11.5664280480244,-4.64172632029319,-2.4340607090218,3.31238894937059,-19.6757175447023,0.841095325921919,0.23905989),
                     c("GO:0051092","positive regulation of NF-kappaB activity",0.862409108843921,-3.43169644603525,5.02299096506065,2.18752072083646,-5.16367588429325,0.926794518496367,0.24020407),
                     c("GO:0042221","response to chemical",22.2704469872048,-4.67343828888502,-4.15994028710556,3.59681693591559,-6.07572071393812,0.902420533275742,0.27670986),
                     c("GO:0001959","regulation of cytokine-mediated signaling pathway",0.839862465475452,3.73670627453505,-5.93970799704764,2.17609125905568,-7.16941133131486,0.832062111484456,0.29253621),
                     c("GO:0060759","regulation of response to cytokine stimulus",0.901865734738741,4.75029706982112,-4.47978815801925,2.20682587603185,-8.29499204066666,0.846334590528637,0.29511386),
                     c("GO:0035457","cellular response to interferon-alpha",0.0563666084211713,-7.30411724043153,0.437006190584048,1.04139268515823,-5.11633856484638,0.875762817568278,0.30802168),
                     c("GO:0097696","receptor signaling pathway via STAT",0.259286398737388,-2.24028125246051,-0.86980437067701,1.67209785793572,-9.80687540164554,0.889798707150165,0.31974088),
                     c("GO:0071360","cellular response to exogenous dsRNA",0.101459895158108,-7.39631323995711,-0.17692696990748,1.27875360095283,-8.97061622231479,0.874996972337628,0.3261111),
                     c("GO:0002831","regulation of response to biotic stimulus",1.94464799053041,4.66407239723672,-5.31212359701012,2.53907609879278,-7.80687540164554,0.832826707915126,0.32611192),
                     c("GO:0043122","regulation of I-kappaB kinase/NF-kappaB signaling",1.40352854968717,3.41372752567685,-5.38029543145439,2.39794000867204,-6.18045606445813,0.828360527978427,0.33636904),
                     c("GO:0019883","antigen processing and presentation of endogenous antigen",0.140916521052928,2.92291737210882,5.35611594159703,1.41497334797082,-14.08512818246,0.777902024372143,0.34575522),
                     c("GO:0032069","regulation of nuclease activity",0.11836987768446,-4.61747303161587,4.67868662917589,1.34242268082221,-6.41793663708829,0.945073159339576,0.37871903),
                     c("GO:0019882","antigen processing and presentation",0.535482780001127,3.6506120401813,4.9663925278114,1.98227123303957,-5.30803489723264,0.862986943089828,0.39727728),
                     c("GO:0080134","regulation of response to stress",7.33329575559439,4.13499510321394,-5.2077883033425,3.11461098423217,-12.2027324591693,0.802985573589072,0.3984021));

one.data <- data.frame(revigo.data);
names(one.data) <- revigo.names;
one.data <- one.data [(one.data$plot_X != "null" & one.data$plot_Y != "null"), ];
one.data$plot_X <- as.numeric( as.character(one.data$plot_X) );
one.data$plot_Y <- as.numeric( as.character(one.data$plot_Y) );
one.data$log_size <- as.numeric( as.character(one.data$log_size) );
one.data$value <- as.numeric( as.character(one.data$value) );
one.data$frequency <- as.numeric( as.character(one.data$frequency) );
one.data$uniqueness <- as.numeric( as.character(one.data$uniqueness) );
one.data$dispensability <- as.numeric( as.character(one.data$dispensability) );


one.data<-dplyr::filter(one.data, dispensability < .25, uniqueness > .844) # change this: filter for highest-confidence/most interesting terms
ggplot(one.data, aes(x = -value, y = reorder(description, -value), 
                     color = log_size)) +
    geom_point(size = 5) +
    theme_classic() +
    scale_color_continuous("GO term\n frequency", limits = c(0,3.6)) +
    labs(y=NULL,x="-log(p)") +
    scale_x_continuous(lim = c(0, 80), breaks = seq(0, 80, 40)) + # change this to fit what scale you need
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.position = "none") # I removed the legend bc it's in another ggplot, but you probably want to include it
ggsave("Revigo_IC3.jpg", 
       width = 7, height = 6.5, 
       units = "in")
