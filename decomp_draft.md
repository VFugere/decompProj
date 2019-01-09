**Title:** Litter decomposition in afrotropical streams: effects of land use, home-field advantages, and terrestrial herbivory

**Authors:** Vincent Fugère, Emily Lostchuck & Lauren J. Chapman

_Target journals: Freshwater Science (IF:3.07), River Research and Application (2.06), Hydrobiologia (2.05), Limnologica (1.8), Inland waters (1.6), Journal of Limnology (1.2), Limnology (1.13), Archiv für hydrobiologie (1.1), African Journal of Aquatic Science (0.67)_

**Abstract**

Land use can strongly affect litter decomposition, a key ecosystem function in low-order streams. Recent evidence suggests that additional drivers of decomposition rates could include 'home-field advantages', when litter decomposes faster at nearby than distant sites, and terrestrial herbivory, whereby inducible plant defenses triggered by herbivore damage can slow down leaf decomposition in streams. To compare the relative importance of these three drivers, we conducted a decomposition experiment in an afrotropical stream system, manipulating land use (farm vs. forest sites), home-field advantage (home vs. away from site of leaf collection), and terrestrial herbivory (using leaves varying in their extent of herbivore damage). We measured decomposition in both fine-mesh and coarse-mesh litter bags to compare drivers of microbial vs. invertebrate-mediated decomposition. Microbial decomposition in fine-mesh bags was unaffected by experimental treatments. In coarse-mesh bags, land use was the only significant (and a strong) driver of decomposition rate, most likely because invertebrate shredders are absent from farm sites. We conclude that home-field advantages and terrestrial herbivory are unimportant drivers of litter decomposition rates in afrotropical streams, at least relative to a major anthropogenic disturbance such as agricultural land use.

**Intro**

Leaf litter breakdown a key ecosystem function contributing x % of carbon cycling in streams. Natural factors driving rates of litter breakdown include x and x. Anthropogenic disturbances affecting one or several of these physico-chemical and biological parameters thus also influence decomposition rates; for example, agricultural land use often slows down litter breakdown considerably (). This response of litter breakdown to disturbance has led some authors to suggest that litter decomposition should be employed as an indicator of stream integrity ().

Recent experimental evidence has revealed hitherto under-appreciated influences on decomposition rates, namely the spatial variation in terrestrial herbivory and in local adaptation of stream communities to role of intraspecific trait variation in trees. These effects have yet to be tested in a variety of stream systems. Moreover, no study has yet quantified the relative importance of these drivers relative to classic natural or anthropogenic drivers.

Here we X. reciprocal transplant experiment. Compare effect sizes.

**Methods**

_Study system_

Fieldwork was conducted in and around Kibale National Park, a 795 km2 mid-altitude (1100-1600 m) rainforest located in southwestern Uganda (0°13' – 0°41' N, 30°19' – 30°32' E). We selected one forested and one agricultural (farm) stream from each of two main watersheds draining the park (the Mpanga River and Dura River watersheds, both of which are subwatersheds of the Nile Basin). All sites are located within 3 km of the Makerere University Biological Field Station in the northwestern part of Kibale. The two study streams inside the park have a fully forested and protected watershed. The two study streams outside of the park have a watershed dominated by intensive agriculture of food and cash crops, pastures for goats and cows, and sparse exotic trees (pine and eucalyptus trees) planted for timber. All sites are small first-order streams (<1.5 m mean wetted width; < 10 cm mean depth) with a similar geomorphology and hydrology, but they vary greatly in water chemistry and community composition based on land use; these differences have been described extensively elsewhere (**FW BIOL**). Briefly, the two farm sites have a much lower canopy cover than the two forest sites, as well as higher water temperature, higher turbidity, lower specific conductance, lower nitrogen and phosphorus concentrations, and a much lower richness and biomass of benthic invertebrates. Invertebrate shredders dominate the composition of forested sites but are largely absent from farm sites (**ecosphere**), leading to slower litter breakdown rates at farm sites.

We delineated a 100-m study reach in each stream to conduct a litterbag experiment (Graça, Bärlocher, & Gessner, 2005). This decomposition experiment focused on the tree species _Neoboutonia macrocalyx_ Pax (Euphorbiaceae), the only species of tree occuring within < 5 m of all four sites. This tree, abundant within the park (**ref**), is kept for shade outside of the park in pastures that are otherwise entirely cleared. N. macrocalyx leaves are readily consumed by both terrestrial herbivores and by macroinvertebrate shredders in streams. Other than _N. macrocalyx_, vegetation in the riparian zone of agricultural sites was composed exclusively of grasses and emergent macrophytes, while the riparian zone of forested sites included a high diversity of trees, shrubs and ferns, with a canopy height generally > 20 m.

_Reciprocal transplant experiment_

We conducted our decomposition experiment in June-August 2011. Twenty leaves varying in their extent of herbivore damage were collected from 3-5 trees around each stream (total = 80 leaves). Leaves were divided into four fragments by cutting two halves along the main stem, and then folding and cutting each half into two additional fragments. This cutting procedure was used to assess whether a given leaf decomposes faster at its home site vs. at a distant site. All fragments (n = 320) were individually-marked, pressed in-between two glass plates, and then photographed against a white background and from a constant distance using a digital camera mounted on a tripod. Leaf fragments were then air-dried for 48 hours in a food dehydrator and weighed individually. Air-dried mass was converted to leaching-adjusted ash-free dry mass (REF) using an air-dried to oven-dried mass conversion equation, a leaching correction factor, and a mean ash content that were previously-calculated for this tree species (REF). Air-dried leaf fragments were randomly assigned to litter bags constructed of 0.5 or 10 mm mesh (fine and coarse-mesh bags, respectively). Decomposition in fine-mesh bags occurs via microbial breakdown only, whereas decomposition in coarse-mesh bags is the product of both microbial and macroinvertebrate-mediated breakdown.

Bags were randomly assigned to arrays of 5 fine-mesh and 5 coarse-mesh bags attached at regular intervals along a 1 m circular metal line. Eight circular arrays were randomly assigned to a study stream, leading to 80 litter bags (40 of each mesh type) being deployed at each site. Arrays were anchored to the stream bottom using stones and twist ties. Every 7 days over 28 days, two arrays (20 bags; 10 of each mesh type) were pulled from each site. The litter content of retrieved bags was rinsed, dried to constant mass for 48 hours at 60 C, weighed, and combusted at 550 C for 4 hours to calculate AFDM of litter remaining in each bag. The proportion of leaf mass that decomposed in the stream, the response variable used in analyses described below, is given by: 1-(AFDM after stream exposure / initial AFDM).

_Terrestrial herbivory_

Leaf photographs were used to quantify damage by terrestrial herbivores, defined as the relative area of a leaf consumed prior to leaf collection. The software Image J version 1.46 (available at https://imagej.nih.gov) was used to measure total and consumed leaf area. An outline of the leaf fragment was first traced assuming no herbivore damage (Fig. 1). This selected area was then binarized using the treshold function in Image J, leading to pixels forming part of the leaf being assigned a brightness value of 1 (white) and pixels contained within holes being assigned a brightness value of 0 (black). The % of original (undamaged) leaf area consumed by terrestrial herbivores is given by: black pixels / (black + white pixels).

Aggregate at leaf and tree level because damage can have both geometrical effects (at fragment scale) and constitutive effects (at tree scale).

_Statistical analyses_

We used R version X for all analyses.

All data and analysis code are available at https://github.com/VFugere/decompProj.

to do: add tree identity to dataset. add random effect of tree. add tree-level damage and test if predicts. Re-run random slopes model with ln (prop) instead of prop as relationship is non-linear.

**Results**

dcsdcs

**Discussion**

dscdcds

**Acknowledgements**

We thank the Uganda Wildlife Authority, the Office of the President of Uganda, and the Uganda National Council of Science and Technology for permission to conduct research in Kibale. The authors are also grateful to Colin Chapman and Jessica Rothman for sharing field and laboratory equipment, to the research assistants of the Kibale Fish and Monkey Project and to the staff at Makarere University Biological Field Station for their help in the field, and to the Natural Sciences and Engineering Research Council of Canada and the Canada Research Chair program for financial support.

**Figure legends**

_Figure 1_. Measurement of herbivore damage. The photograph shows a N. macrocalyx leaf fragment with some herbivore damage. The solid black line depicts leaf area assuming no damage while dashed lines encircle areas eaten by herbivores. Leaf damage was quantified as the total number of pixels inside dashed polygons divided by the number of pixels within the solid polygon.

_Figure 2_. Decomposition rates in fine-mesh (a) and coarse-mesh (b) bags as a function of land use (left panels), home-field advantage (middle panels), and terrestrial herbivory (right panels). Lines and shaded polygons indicate means +/- 95% confidence intervals of the mean. For terrestrial herbivory, data are shown for the final time point of the experiment.

_Figure 3_. Results of GLMMs quantifying the influence of five variables on the proportion of leaf litter decomposed by the end of the experiment. Separate models were fitted for fine-mesh bags (a) and coarse-mesh bags (b). Symbols and error bars indicate parameter estimates +/- 95 % confidence intervals; statistically-significant effects with confidence intervals that do not overlap zero are shown in solid colour, while non-significant effects are shown in transparent colour. Models were fitted using data from the final time point of the experiment.
