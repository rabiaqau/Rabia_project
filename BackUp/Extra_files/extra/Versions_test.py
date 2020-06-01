

from ROOT import *
import AtlasStyle
import AtlasUtils
import ROOT
import math
import numpy as np
import os.path
import pickle
import os.path
import pickle
import sys
import glob
from functions import *





run_number= [ 348354,  348354,348403,348495,348495,348511,348511,348534,
  348534,348534,348534,348534,348534,348610,348618,349011,
  349011,349014,349033,349033,349051,349051,349111,349114,
  349114,349169,349169,349268,349268,349309,349327,349335,
  349451,349481,349481,349533,349592,349637,349637,349646,
  349693,349693,349841,349842,349944,349944,349977,350013,
  350013,350067,350121,350121,350144,350144,350160,350184,
  350184,350220,350310,350361,350361,350440,350440,350479,
  350479,350531,350531,350676,350676,350676,350676,350682,
  350682,350842,350880,350880,351062,351223,351223,351364,
  351364,351550,351628,351671,351671,351969,351969,352107,
  352274,352448,354124,354359,354396,355261,355261,355331,
  355389,355416,355446,355446,355468,355468,355529,355544,
  355754,355754,355995,355995,356177,356250,357193,357283,
  357409,357409,357451,357500,357539,357620,357620,357679,
  357713,357713,357772,357821,357887,357887,358031,358031,
  358215,358300,358395,358395,358516,358541,358541,358615,
  358615,358656,358985,359058,359124,359170,359171,359191,
  359191,359286,359310,359355,359472,359541,359541,359678,
  359717,359735,359735,359823,359823,359872,359918,360026,
  360063,360161,360209,360244,360309,360309,360402,360414,
  361738,361738,361795,361862,361862,362354,362445,362445,
  362552,362661,362661,363033,363129,363198,363198,363400,
  363664,363710,363710,363738,363830,363830,363910,363947,
  363979,364030,364030,364076,364098,364098,364214,364292,
  364292]








print len(run_number),"ru"
fom_obtained_old_ver_trk= [ 0.94617604 ,0.99467061 ,0.99271133 ,0.9634011 , 0.98808589 ,0.9823237
 ,0.9910357 , 1.07930081 ,0.97179523 ,1.00898122 ,0.98976677 ,1.04782714
 ,1.0555301 , 0.98427902 ,0.97476947 ,0.97986116 ,0.99302172 ,0.99744707
 ,0.97448691 ,0.99582535 ,0.98261557 ,1.00487562 ,0.97738088 ,0.97498187
 ,0.97578745 ,0.96516227 ,0.98662094 ,0.98243277 ,0.98399371 ,0.98629504
 ,0.97854685 ,0.97835804 ,0.98935786 ,0.97928045 ,1.00045768 ,0.97268002
 ,0.9894229 , 0.98893306 ,0.99599781 ,0.97035258 ,0.98042202 ,1.01358508
 ,0.97096185 ,0.99189686 ,0.97435702 ,1.00006068 ,0.98552604 ,0.97984397
 ,0.99577324 ,0.98010426 ,0.9867755 , 1.01480941 ,0.97436261 ,0.98143808
 ,0.974724  ,0.9797575 , 0.98340875 ,0.97720321 ,0.97874085 ,0.98375181
 ,1.00386421 ,0.98511142 ,1.00098836 ,0.96333882 ,1.00134925 ,0.9615018
 ,0.99849223 ,1.09486305 ,1.07814984 ,1.10020445 ,1.08007596 ,0.97267046
 ,0.99251164 ,0.96901926 ,0.97640376 ,1.03234203 ,0.96759016 ,0.97940287
 ,0.98129059 ,0.96784573 ,0.99620349 ,0.96709034 ,0.98624317 ,0.96937256
 ,0.99911473 ,0.97655927 ,0.98679374 ,0.98176761 ,0.98212432 ,0.97714577
 ,0.99201736 ,0.96843153 ,0.98810756 ,0.97044058 ,0.98861871 ,0.99881801
 ,0.97429991 ,0.97610346 ,1.09912676 ,1.09999281 ,0.96688769 ,1.00659493
 ,0.96072795 ,0.96885136 ,0.98191897 ,1.00181527 ,0.97688621 ,0.99580785
 ,0.96806665 ,0.97475231 ,0.99618464 ,0.98385831 ,0.96009425 ,0.9903136
 ,0.99146625 ,0.97835837 ,0.97888749 ,0.97655343 ,0.97073458 ,0.98927626
 ,0.97102884 ,0.97966079 ,0.9699273 , 0.98869689 ,0.94847457 ,0.96608551
 ,0.97667841 ,0.98573271 ,0.9526237 , 0.96212437 ,0.96387783 ,1.01775155
 ,0.98831575 ,0.98169125 ,0.98371529 ,0.97148114 ,0.99580995 ,0.99143252
 ,0.97940126 ,0.94311887 ,0.99492833 ,0.95910967 ,0.99503033 ,0.93932317
 ,0.99871086 ,0.98628705 ,0.99892248 ,0.96666241 ,1.00404753 ,0.95754567
 ,0.99233235 ,0.96626615 ,0.99384098 ,0.96786871 ,1.02281571 ,0.97331537
 ,0.98292205 ,0.99555973 ,0.97321595 ,0.98546494 ,0.96047788 ,0.96563581
 ,0.96842561 ,0.99169914 ,0.97482011 ,0.97936566 ,0.95896447 ,0.96339782
 ,0.96483982 ,0.97859508 ,0.97724133 ,0.97085679 ,0.99049099 ,0.96598038
 ,0.9696226 , 0.97240673 ,0.98151256 ,0.96889311 ,0.9929445 , 1.00095763
 ,0.96796748 ,0.96432741 ,0.99379807 ,0.98219361 ,0.98116615 ,0.97218842
 ,0.99140833 ,0.99142439 ,0.9676543 , 1.00295925 ,0.9737063 , 0.97581794
 ,1.00990077 ,0.96495414 ,0.98832661 ,0.98375239 ,0.97326062 ,0.98217985
 ,1.01062475 ,0.97495042 ,1.04984444]


print len(fom_obtained_old_ver_trk),"old track"


new_ver_trk_fom= [ 0.98408144 ,0.94617604 ,0.9823237 , 1.07930081 ,0.97179523 ,1.00898122
 ,0.97986116 ,0.97448691 ,0.98261557 ,0.97738088 ,0.97498187 ,0.96516227
 ,0.98243277 ,0.98629504 ,0.97854685 ,0.97835804 ,0.98935786 ,0.97928045
 ,0.97268002 ,0.9894229 , 0.98893306 ,0.97035258 ,0.98042202 ,0.97096185
 ,0.97435702 ,0.98552604 ,0.97984397 ,0.98010426 ,0.9867755 , 0.97436261
 ,0.9797575 , 0.97720321 ,0.98375181 ,0.98511142 ,0.96333882 ,0.9615018
 ,1.09486305 ,0.97267046 ,0.96901926 ,0.97640376 ,0.96759016 ,0.97940287
 ,0.96784573 ,0.96709034 ,0.96937256 ,0.97655927 ,0.98176761 ,0.98212432
 ,0.96843153 ,0.97044058 ,0.97429991 ,0.97610346 ,1.09912676 ,1.09999281
 ,0.96688769 ,0.96072795 ,0.96885136 ,0.98191897 ,0.97688621 ,0.96806665
 ,0.97475231 ,0.98385831 ,0.96009425 ,0.97835837 ,0.97655343 ,0.97102884
 ,0.94847457 ,0.97667841 ,0.9526237 , 0.96212437 ,0.96387783 ,0.98169125
 ,0.97148114 ,0.97940126 ,0.94311887 ,0.95910967 ,0.93932317 ,0.98628705
 ,0.96666241 ,0.95754567 ,0.96626615 ,0.96786871 ,0.97331537 ,0.97321595
 ,0.96047788 ,0.96563581 ,0.96842561 ,0.97482011 ,0.95896447 ,0.96483982
 ,0.97085679 ,0.96598038 ,0.9696226 , 0.98151256 ,0.96889311 ,0.96432741
 ,0.97218842 ,0.9676543 , 0.97581794 ,0.96495414 ,0.97326062 ,0.97495042
 ,1.04747521 ,0.99467061 ,0.99271133 ,0.9634011 , 0.98808589 ,0.9910357
 ,0.98976677 ,1.04782714 ,1.0555301 , 0.98427902 ,0.97476947 ,0.99302172
 ,0.99744707 ,0.99582535 ,1.00487562 ,0.97578745 ,0.98662094 ,0.98399371
 ,1.00045768 ,0.99599781 ,1.01358508 ,0.99189686 ,1.00006068 ,0.99577324
 ,1.01480941 ,0.98143808 ,0.974724  ,0.98340875 ,0.97874085 ,1.00386421
 ,1.00098836 ,1.00134925 ,0.99849223 ,1.07814984 ,1.10020445 ,1.08007596
 ,0.99251164 ,1.03234203 ,0.98129059 ,0.99620349 ,0.98624317 ,0.99911473
 ,0.98679374 ,0.97714577 ,0.99201736 ,0.98810756 ,0.98861871 ,0.99881801
 ,1.00659493 ,1.00181527 ,0.99580785 ,0.99618464 ,0.9903136 , 0.99146625
 ,0.97888749 ,0.97073458 ,0.98927626 ,0.97966079 ,0.9699273 , 0.98869689
 ,0.96608551 ,0.98573271 ,1.01775155 ,0.98831575 ,0.98371529 ,0.99580995
 ,0.99143252 ,0.99492833 ,0.99503033 ,0.99871086 ,0.99892248 ,1.00404753
 ,0.99233235 ,0.99384098 ,1.02281571 ,0.98292205 ,0.99555973 ,0.98546494
 ,0.99169914 ,0.97936566 ,0.96339782 ,0.97859508 ,0.97724133 ,0.99049099
 ,0.97240673 ,0.9929445 , 1.00095763 ,0.96796748 ,0.99379807 ,0.98219361
 ,0.98116615 ,0.99140833 ,0.99142439 ,1.00295925 ,0.9737063 , 1.00990077
 ,0.98832661 ,0.98375239 ,0.98217985 ,1.01062475 ,1.04984444]
print len(new_ver_trk_fom)

































# fom_obtained_from_old_ver_lcd= [ 0.94479077 ,0.96981804 ,0.98310484 ,0.96960838 ,0.99721309 ,0.97859678
#  ,0.98277844 ,0.98353102 ,0.97295216 ,0.97779819 ,0.99053743 ,0.9941207
#  ,0.98159251 ,0.98157162 ,0.96410473 ,0.97740462 ,0.99024774 ,0.98308931
#  ,0.9703834 , 0.99218616 ,0.98272394 ,0.98106699 ,0.97467586 ,0.9646292
#  ,0.96118025 ,0.96265077 ,0.98201892 ,0.97970973 ,0.97437146 ,0.9832676
#  ,0.96771365 ,0.96630521 ,0.98939387 ,0.97678237 ,0.98109476 ,0.92594079
#  ,0.98896493 ,0.98999214 ,0.99119853 ,0.93358087 ,0.98337649 ,1.00536013
#  ,0.97327249 ,0.99556396 ,0.99614453 ,0.99286389 ,0.97868863 ,0.97516298
#  ,0.99266974 ,0.97769896 ,0.98201399 ,0.98631922 ,0.96939471 ,0.97632729
#  ,0.99436839 ,0.97348945 ,0.9751248 , 0.96850708 ,0.97353115 ,0.97972777
#  ,1.00197001 ,0.98199297 ,1.00717436 ,0.9611985 , 0.99547831 ,0.95768551
#  ,0.97795445 ,1.09979085 ,1.05567348 ,1.11556023 ,1.09491647 ,0.97190854
#  ,0.99708213 ,0.96664158 ,0.96738488 ,1.03158145 ,0.96681293 ,0.97740633
#  ,0.98041496 ,0.9690278 , 1.0009411 , 0.96526154 ,0.98430365 ,0.96520706
#  ,0.98448012 ,0.97502548 ,0.98467486 ,0.97661548 ,0.98098743 ,0.98472257
#  ,0.9670974 , 0.96002138 ,0.98335001 ,0.96384512 ,0.98394532 ,0.99691019
#  ,0.97109458 ,0.97484678 ,1.08446973 ,1.08745203 ,0.9613122 , 1.00998176
#  ,0.95792207 ,0.96615245 ,0.97848424 ,0.99489338 ,0.97443147 ,0.9986848
#  ,0.96495711 ,0.93399229 ,0.99740485 ,0.98236511 ,0.96002349 ,0.98992306
#  ,0.98497221 ,0.97166061 ,0.97758743 ,0.97689411 ,0.97333211 ,0.99771486
#  ,0.97929887 ,0.98169094 ,0.96967752 ,0.99238615 ,0.95102362 ,0.96960969
#  ,0.97531257 ,0.98957207 ,0.95123119 ,0.95524386 ,0.95869742 ,1.0132884
#  ,0.99023562 ,0.97355362 ,0.98122455 ,0.96583331 ,0.99503277 ,0.99207951
#  ,0.96767391 ,0.91314687 ,0.98864261 ,0.95136336 ,0.98937526 ,0.93489532
#  ,1.00041932 ,0.98354584 ,0.99829888 ,0.96114613 ,1.00315652 ,0.95319634
#  ,0.99464404 ,0.96433482 ,0.99469847 ,0.92497534 ,1.00671062 ,0.97299562
#  ,0.98151828 ,1.0200453 , 0.96718632 ,0.9845276 , 0.95620283 ,0.96054905
#  ,0.96832439 ,0.99271968 ,0.97320466 ,0.97902808 ,0.95406077 ,0.99201298
#  ,0.9608381 , 0.97847994 ,0.97230516 ,0.96575292 ,0.99158029 ,0.95973068
#  ,0.96621662 ,0.97357031 ,0.98021159 ,0.96337436 ,0.98927464 ,1.00359326
#  ,0.96598458 ,0.96271978 ,0.98936624 ,0.97935158 ,0.97686167 ,0.94235194
#  ,0.99011515 ,0.99193613 ,0.96239582 ,0.99898778 ,0.97120157 ,0.97172501
#  ,0.99813389 ,0.9620458 , 0.98732274 ,0.98201493 ,0.97083844 ,0.98443647
#  ,1.012307 , 0.97054109 ,0.99851569]


# new_ver_lcd_fom= [ 1.03565269 ,0.94479077 ,0.97859678 ,0.98353102 ,0.97295216 ,0.97779819
#  ,0.97740462 ,0.9703834 , 0.98272394 ,0.97467586 ,0.9646292 , 0.96265077
#  ,0.97970973 ,0.9832676 , 0.96771365 ,0.96630521 ,0.98939387 ,0.97678237
#  ,0.92594079 ,0.98896493 ,0.98999214 ,0.93358087 ,0.98337649 ,0.97327249
#  ,0.99614453 ,0.97868863 ,0.97516298 ,0.97769896 ,0.98201399 ,0.96939471
#  ,0.97348945 ,0.96850708 ,0.97972777 ,0.98199297 ,0.9611985 , 0.95768551
#  ,1.09979085 ,0.97190854 ,0.96664158 ,0.96738488 ,0.96681293 ,0.97740633
#  ,0.9690278 , 0.96526154 ,0.96520706 ,0.97502548 ,0.97661548 ,0.98098743
#  ,0.96002138 ,0.96384512 ,0.97109458 ,0.97484678 ,1.08446973 ,1.08745203
#  ,0.9613122 , 0.95792207 ,0.96615245 ,0.97848424 ,0.97443147 ,0.96495711
#  ,0.93399229 ,0.98236511 ,0.96002349 ,0.97166061 ,0.97689411 ,0.97929887
#  ,0.95102362 ,0.97531257 ,0.95123119 ,0.95524386 ,0.95869742 ,0.97355362
#  ,0.96583331 ,0.96767391 ,0.91314687 ,0.95136336 ,0.93489532 ,0.98354584
#  ,0.96114613 ,0.95319634 ,0.96433482 ,0.92497534 ,0.97299562 ,0.96718632
#  ,0.95620283 ,0.96054905 ,0.96832439 ,0.97320466 ,0.95406077 ,0.9608381
#  ,0.96575292 ,0.95973068 ,0.96621662 ,0.98021159 ,0.96337436 ,0.96271978
#  ,0.94235194 ,0.96239582 ,0.97172501 ,0.9620458 , 0.97083844 ,0.97054109
#  ,1.01250532 ,0.96981804 ,0.98310484 ,0.96960838 ,0.99721309 ,0.98277844
#  ,0.99053743 ,0.9941207 , 0.98159251 ,0.98157162 ,0.96410473 ,0.99024774
#  ,0.98308931 ,0.99218616 ,0.98106699 ,0.96118025 ,0.98201892 ,0.97437146
#  ,0.98109476 ,0.99119853 ,1.00536013 ,0.99556396 ,0.99286389 ,0.99266974
#  ,0.98631922 ,0.97632729 ,0.99436839 ,0.9751248 , 0.97353115 ,1.00197001
#  ,1.00717436 ,0.99547831 ,0.97795445 ,1.05567348 ,1.11556023 ,1.09491647
#  ,0.99708213 ,1.03158145 ,0.98041496 ,1.0009411 , 0.98430365 ,0.98448012
#  ,0.98467486 ,0.98472257 ,0.9670974 , 0.98335001 ,0.98394532 ,0.99691019
#  ,1.00998176 ,0.99489338 ,0.9986848 , 0.99740485 ,0.98992306 ,0.98497221
#  ,0.97758743 ,0.97333211 ,0.99771486 ,0.98169094 ,0.96967752 ,0.99238615
#  ,0.96960969 ,0.98957207 ,1.0132884 , 0.99023562 ,0.98122455 ,0.99503277
#  ,0.99207951 ,0.98864261 ,0.98937526 ,1.00041932 ,0.99829888 ,1.00315652
#  ,0.99464404 ,0.99469847 ,1.00671062 ,0.98151828 ,1.0200453 , 0.9845276
#  ,0.99271968 ,0.97902808 ,0.99201298 ,0.97847994 ,0.97230516 ,0.99158029
#  ,0.97357031 ,0.98927464 ,1.00359326 ,0.96598458 ,0.98936624 ,0.97935158
#  ,0.97686167 ,0.99011515 ,0.99193613 ,0.99898778 ,0.97120157 ,0.99813389
#  ,0.98732274 ,0.98201493 ,0.98443647 ,1.012307  ,0.99851569]
# print len(new_ver_lcd_fom)
