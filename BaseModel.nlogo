;; THIS IS MAIN FILE. MOST procedure called from here
__includes
[
  ;# PATH 4.0
  "PATH4main.nls";## Main

  "global-initialization.nls" "main-code-helper-functions.nls" ;helper functions;## GLOBAL INITIALIZATIONS
  "initial-people-ECNA.nls" "breed-susceptibles.nls" "dryrun.nls" ;;## AGENT / GLOBOCAL INITALIZATIONS

  "disease_progression.nls" "Data.nls" "testing-frequency.nls" "set-dropout-care.nls" "manage-care-continuum_juri.nls" ;; ## HIV disease progression
  "breed-T-tree-links.nls" "breed-people.nls"

  "transmissions-ECNA.nls";## Transmisison -  Bernoulli
  "ECNA-generateNetwork-v2.nls" "contact-activation.nls" ;## ECNA ALGORITHM
  "age-functions.nls" "jurisdiction.nls" "riskGroupPreference.nls" ; ## DEMOGRAPHICS risk groups selections

  ;# WRITE_OUPU0TS
  "write-mixing.nls" "write-epidemic-features-UI.nls" ; PATH 4.0
  "write-output.nls" "write-output-headers.nls" "Visualization.nls"  ; PATH 2.0 and 3.0 (useful for sensitivity analysis setup)

  "overrides_example.nls" "nooverride.nls" ;needed for PATH UI - output visualization app


  ;# Multiple diseases
  "write-HPV.nls" "multiple-disease.nls" "debug-multi-disease.nls"

]

;; Extensions of NetLogo used
extensions [py py matrix csv profiler table rnd]


;;PATH4.0
to runExperiment
  override;; needed for PATHUI- output visualization app
  nooverride ;; needed for PATHUI- output visualization app
  reset-timer
  random-seed 1076
  let run_num 1
  while [run_num <= maxRun] [
    setupECNA
    repeat termination-ticks[runECNA  if count people with [min nodes-link-active-age <= 13] > 0  [print "age below"] ]
    print (list "completed run" run_num "of"  maxRun "; timer" timer)
    set run_num run_num + 1
  ]

end
@#$#@#$#@
GRAPHICS-WINDOW
198
13
805
474
-1
-1
10.512
1
10
1
1
1
0
0
0
1
-28
28
-21
21
1
1
1
ticks
30.0

PLOT
439
811
756
949
CD4 count at diagnosis - Heterosexuals
Time
Proportion
0.0
10.0
0.0
1.0
true
true
"" ""
PENS
"<=200" 1.0 0 -16777216 true "" "plot count people with [CD4-diagnosis <= 200 and stage > 2 and sex <= 2 and dead = 0] / count people with [stage > 2 and sex <= 2 and dead = 0]"
"200-350" 1.0 0 -13345367 true "" "plot count people with [CD4-diagnosis <= 350 and CD4-diagnosis > 200 and stage > 2 and sex <= 2 and dead = 0] / count people with [stage > 2 and sex <= 2 and dead = 0]"
"350-500" 1.0 0 -5825686 true "" "plot count people with [CD4-diagnosis <= 500 and CD4-diagnosis > 350 and stage > 2 and sex <= 2 and dead = 0] / count people with [stage > 2 and sex <= 2 and dead = 0]"
">500" 1.0 0 -10899396 true "" "plot count people with [CD4-diagnosis > 500 and stage > 2 and sex <= 2 and dead = 0] / count people with [stage > 2 and sex <= 2 and dead = 0]"

INPUTBOX
1037
101
1099
161
time-unit
12.0
1
0
Number

PLOT
0
812
440
1070
Distribution of people living with HIV/AIDS
Time
Proportion
0.0
10.0
0.0
1.0
true
true
"" ""
PENS
"Heterosexual-female" 1.0 0 -13345367 true "" "plot count people with [sex = 1 and dead = 0]"
"MSM" 1.0 0 -13791810 true "" "plot count people with [sex = 3 and dead = 0]"
"Heterosexual-male" 1.0 0 -10899396 true "" "plot count people with [sex = 2 and dead = 0]"

PLOT
757
947
1132
1073
Stage Distribution- MSM
Time
Proportion
0.0
10.0
0.0
1.0
true
true
"" ""
PENS
"acute-unaware" 1.0 0 -16777216 true "" "plot count people with [stage = 1 and dead = 0 and sex = 3] / count people with [dead = 0 and sex = 3]"
"non-acute-unaware" 1.0 0 -7500403 true "" "plot count people with [stage = 2 and dead = 0 and sex = 3] / count people with [dead = 0 and sex = 3]"
"aware-no care" 1.0 0 -2674135 true "" "plot count people with [stage = 3 and dead = 0 and sex = 3] / count people with [dead = 0 and sex = 3]"
"aware-care-no ART" 1.0 0 -955883 true "" "plot count people with [stage = 4 and dead = 0 and sex = 3] / count people with [dead = 0 and sex = 3]"
"ART-not suppressed" 1.0 0 -1184463 true "" "plot count people with [stage = 5 and dead = 0 and sex = 3] / count people with [dead = 0 and sex = 3]"
"ART-suppressed" 1.0 0 -10899396 true "" "plot count people with [stage = 6 and dead = 0 and sex = 3] / count people with [dead = 0 and sex = 3]"

PLOT
440
948
758
1068
CD4 count at diagnosis - MSM
Time
Proportion
0.0
10.0
0.0
1.0
true
true
"" ""
PENS
"<=200" 1.0 0 -16777216 true "" "plot count people with [CD4-diagnosis <= 200 and stage > 2 and sex = 3 and dead = 0] / count people with [stage > 2 and sex = 3 and dead = 0]"
"200-350" 1.0 0 -13345367 true "" "plot count people with [CD4-diagnosis <= 350 and CD4-diagnosis > 200 and stage > 2 and sex = 3 and dead = 0] / count people with [stage > 2 and sex = 3 and dead = 0]"
"350-500" 1.0 0 -5825686 true "" "plot count people with [CD4-diagnosis <= 500 and CD4-diagnosis > 350 and stage > 2 and sex = 3 and dead = 0] / count people with [stage > 2 and sex = 3 and dead = 0]"
">500" 1.0 0 -10899396 true "" "plot count people with [CD4-diagnosis > 500 and stage > 2 and sex = 3 and dead = 0] / count people with [stage > 2 and sex = 3 and dead = 0]"

MONITOR
1222
1027
1272
1072
undiag
count people with [stage = 2 and dead = 0 and index? = true and new-infections = 0]
17
1
11

INPUTBOX
1846
799
1896
859
goal
1.0
1
0
Number

INPUTBOX
771
653
869
713
maxDegree
128.0
1
0
Number

MONITOR
1138
934
1266
979
NIL
eligible-agents-count
0
1
11

MONITOR
1138
981
1268
1026
NIL
eligible-nonagents-count
0
1
11

PLOT
1137
796
1474
932
Proportion of times agents selected
NIL
NIL
0.0
360.0
0.0
1.0
false
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot  selected-agents / (selected-agents + selected-nonagents)"

MONITOR
1142
1028
1224
1073
NIL
trackTotalPop
17
1
11

CHOOSER
928
115
1029
160
first-year
first-year
2006 2010
0

INPUTBOX
826
176
921
236
simulation-years
21.0
1
0
Number

INPUTBOX
825
109
923
169
initial-infected
500.0
1
0
Number

INPUTBOX
929
176
1029
236
termination-ticks
360.0
1
0
Number

TEXTBOX
1208
771
1396
799
ADHOC MONITORING-PATH4.0
11
0.0
1

TEXTBOX
1802
774
1952
792
PATH 2.0 and PATH 3.0 Inputs
11
0.0
1

TEXTBOX
836
85
986
104
PATH 4.0 INPUTS
15
0.0
1

TEXTBOX
1157
19
1592
403
PATH 4.0 input \n*simulation-years: number of years to simulate from year 2006;\n*maxRun: number of simulation iterations\n*initial-infected: HIV population size prior to dryrun. \n*dry_run_1_duration: dry run for intiating network dynamics, clock will not update during this period\n(duration of dry-run 2: defined in model as = termination ticks - (simulation-years * time-unit) dry run for initializating epidemic and network feature dynamics; clock will update during this period.) \n*time-unit: 12 implies time step is monthly\n*first-year: first year of simulation after both dry runs.\n\nDefaults: \n*first-year: 2006\n*simulation-years: 11 (so will simulate 2006 to 2017)\n*termination-ticks: 240\n*dry_run_1_duration: 50\n*initial-infected: 1000 (take ~30 mins for 30 runs parallel run using Netlogo BehaviorSPace on Intel(R) Core(TM) i9-10900X CPU @ 3.70GHz   3.70 GHz).(Generates ~3500 HIV infected persons by end of simulation). For low computational testing, set value to small number ~200, although sometimes this can create an error if sample is too small for intended distribution, however, probbaility of error is low for 200, so if error, just re-run.
13
0.0
1

INPUTBOX
930
241
1099
301
dry_run_1_duration
50.0
1
0
Number

TEXTBOX
956
46
1106
64
<---- RUN PATH 4.0
13
0.0
1

INPUTBOX
698
719
795
779
num-jurisdictions
1.0
1
0
Number

INPUTBOX
805
718
910
778
juri_model_num
1.0
1
0
Number

SWITCH
599
680
760
713
concurrency-model
concurrency-model
0
1
-1000

BUTTON
822
41
934
74
NIL
runExperiment
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

INPUTBOX
1619
129
1721
189
temp-time-unit
2.0
1
0
Number

INPUTBOX
1620
195
1719
255
increased-risk
1.0
1
0
Number

BUTTON
1604
46
1819
79
NIL
runExperiment-multi-disease
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

INPUTBOX
1727
193
1838
253
model-screening
1.0
1
0
Number

INPUTBOX
1840
194
1951
254
model-vaccine
1.0
1
0
Number

TEXTBOX
610
625
866
663
Network features. Keep at default
15
0.0
1

INPUTBOX
1036
175
1101
235
maxRun
1.0
1
0
Number

INPUTBOX
581
720
694
780
jurisdiction-mixing-within
0.77
1
0
Number

TEXTBOX
1838
52
2039
86
<---- RUN PATH 4.0-MULTI-DISEASE\n
14
0.0
1

TEXTBOX
1659
98
1965
138
PATH 4.0-MULTI-DISEASE INPUTS
16
0.0
1

PLOT
760
805
1131
945
Stage Distribution - Heterosexuals
Time
Proportion
0.0
10.0
0.0
1.0
true
true
"" ""
PENS
"acute-unaware" 1.0 0 -16777216 true "" "plot count people with [stage = 1 and dead = 0 and sex <= 2] / count people with [dead = 0 and sex <= 2]"
"non-acute-unaware" 1.0 0 -7500403 true "" "plot count people with [stage = 2 and dead = 0 and sex <= 2] / count people with [dead = 0 and sex <= 2]"
"aware-no-care" 1.0 0 -2674135 true "" "plot count people with [stage = 3 and dead = 0 and sex <= 2] / count people with [dead = 0 and sex <= 2]"
"aware-care-no ART" 1.0 0 -955883 true "" "plot count people with [stage = 4 and dead = 0 and sex <= 2] / count people with [dead = 0 and sex <= 2]"
"ART-not suppressed" 1.0 0 -1184463 true "" "plot count people with [stage = 5 and dead = 0 and sex <= 2] / count people with [dead = 0 and sex <= 2]"
"ART-suppressed" 1.0 0 -10899396 true "" "plot count people with [stage = 6 and dead = 0 and sex <= 2] / count people with [dead = 0 and sex <= 2]"

@#$#@#$#@
##MOLECULAR PROJECT
1. create-Ttree-link-with in initialize-transmissions.nls



## WHAT IS IT?

This section could give a general understanding of what the model is trying to show or explain.

## HOW IT WORKS

This section could explain what rules the agents use to create the overall behavior of the model.

## HOW TO USE IT

This section could explain how to use the model, including a description of each of the items in the interface tab.  
Initially hit "setup" button on the interface tab

1. It sets up global variables by calling "setup-globals"
	a. globals are all the data that remains the same for entire population

2. It calls "setup-people"
 	a. It generates number of people as mentioned in the "number-people" input button. Sets x% (currently 100%) of the initial people as index-patients. Which means everyone gets infected in quarter 1, i.e., tick 1 (but are un-infected in tick 0 for ease of computation). 

 	b. Initializes all the variable used in the model by calling procedure "set-infected-variables"

Next hit the "go" button. This begins the simulation

1. For "simulation-years" updates all the variables of (dead =0, i.e, non-dead) people in the simulation by calling "update-simulation"

2. Each quarter, i.e., each tick, creates new infections (currently creates new people because 100% of initial population was infected in the beginning and can be changed as needed). Number of new infections = total number of transmissions generated by all the index patients only. 
	a. This is done by calling the procedure "setup-new-transmissions [number]".

	b. The procedure sets the new infections as index-patients = false and calls "set-infected-variables" to initialize the variables. From the next quarter they are included in the simulation, that is, "update-simulation" automatically applies to everyone in the simulation

3. when a person dies, variable "dead = 1" but we do not assign "die" which is in-built netlogo function for permanently removing person from the simulation. This is done so that the population statistics can be collected in the end (using die erases persons data)

After simulation stops hit "display-values"  
1. Displays the required data for index and transmissions separately.  
2. For each peron, at time of death, "set-life-variables" which keeps track of life variables. When display-values is clicked the simulation is modeled to display some of these values. Note: Currently all life-variables are for index-patients ONLY. Change as needed in procedures "set-dead", "set-life-variables", and "display-values"

## THINGS TO NOTICE

This section could give some ideas of things for the user to notice while running the model.

## THINGS TO TRY

This section could give some ideas of things for the user to try to do (move sliders, switches, etc.) with the model.

## EXTENDING THE MODEL

This section could give some ideas of things to add or change in the procedures tab to make the model more complicated, detailed, accurate, etc.

## NETLOGO FEATURES

This section could point out any especially interesting or unusual features of NetLogo that the model makes use of, particularly in the Procedures tab.  It might also point out places where workarounds were needed because of missing features.

## RELATED MODELS

This section could give the names of models in the NetLogo Models Library or elsewhere which are of related interest.

## CREDITS AND REFERENCES

This section could contain a reference to the model's URL on the web if it has one, as well as any other necessary credits or references.
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

die 1
false
0
Rectangle -7500403 true true 45 45 255 255
Circle -16777216 true false 129 129 42

die 2
false
0
Rectangle -7500403 true true 45 45 255 255
Circle -16777216 true false 69 69 42
Circle -16777216 true false 189 189 42

die 3
false
0
Rectangle -7500403 true true 45 45 255 255
Circle -16777216 true false 69 69 42
Circle -16777216 true false 129 129 42
Circle -16777216 true false 189 189 42

die 4
false
0
Rectangle -7500403 true true 45 45 255 255
Circle -16777216 true false 69 69 42
Circle -16777216 true false 69 189 42
Circle -16777216 true false 189 69 42
Circle -16777216 true false 189 189 42

die 5
false
0
Rectangle -7500403 true true 45 45 255 255
Circle -16777216 true false 69 69 42
Circle -16777216 true false 129 129 42
Circle -16777216 true false 69 189 42
Circle -16777216 true false 189 69 42
Circle -16777216 true false 189 189 42

die 6
false
0
Rectangle -7500403 true true 45 45 255 255
Circle -16777216 true false 84 69 42
Circle -16777216 true false 84 129 42
Circle -16777216 true false 84 189 42
Circle -16777216 true false 174 69 42
Circle -16777216 true false 174 129 42
Circle -16777216 true false 174 189 42

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
0
Rectangle -7500403 true true 151 225 180 285
Rectangle -7500403 true true 47 225 75 285
Rectangle -7500403 true true 15 75 210 225
Circle -7500403 true true 135 75 150
Circle -16777216 true false 165 76 116

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.2.2
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="test-experiment" repetitions="1" runMetricsEveryStep="true">
    <setup>setupECNA</setup>
    <go>runMultiDisease</go>
    <final>write-multi-disease 1</final>
    <enumeratedValueSet variable="initial-infected">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="time-unit">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="temp-time-unit">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="termination-ticks">
      <value value="360"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="simulation-years">
      <value value="21"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model-screening">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model-vaccine">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="increased-risk">
      <value value="1"/>
    </enumeratedValueSet>
    <steppedValueSet variable="random_seed" first="1" step="1000" last="10000"/>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
