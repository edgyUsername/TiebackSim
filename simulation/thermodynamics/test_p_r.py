import setup_comp_data as setPVT
P=2e5
T=-273.15+300
##pvt={'nitrogen': {'comp': 0.79, 'params': {'acc': 0.0377, 'mm': 0.028, 't_c': 126.2, 'p_c': 34, 'v_c': 0.08921, '_id': 'nitrogen'}}, 'oxygen': {'comp': 0.21, 'params': {'acc': 0.0222, 'mm': 0.032, 't_c': 154.58, 'p_c': 50.43, 'v_c': 0.0734, '_id': 'oxygen'}}}

pvt= setPVT.get_props_and_equalize([(63,'octane'),(36,'methane')])

######################fugacity calc from PR EoS########################
##print pr.get_wilson_K(pvt,T,P)
##print pr.find_vapor_frcn(pvt,T,P,pr.get_wilson_K(pvt,T,P))

######################cubic solver test########################
##print '1 real:\t',pr.cubic_solver(1,1,1,-3)
##print '3 real, 3 the same:\t',pr.cubic_solver(1,-3,3,-1)
##print '3 real, 2 the same:\t',pr.cubic_solver(1,-5,8,-4)
##print '3 real, distinct:\t',pr.cubic_solver(1,-6,11,-6)

######################mixing rules test########################
##print pr.van_der_waal_coefs(pvt,T+273.14,P,comp=None)

######################fugacity calc from PR EoS########################
##import math
##pvt_fug_test=setPVT.get_props_and_equalize([(2,'ethylene'),(8,'argon')])
##P_fug_test=125e5
##T_fug_test=24.95
##lnfug,V= pr.fug_minimum_gibbs(pvt_fug_test,T_fug_test,P_fug_test,1,returnV=True)
##print "P= %s atm, T= %s C" %(P_fug_test/1.013e5,T_fug_test)
##print V*100*100*100
##for i in lnfug:
##    print i,'\t',math.exp(lnfug[i])

##################### lnK calc from EoS#############################
##import math
##lnK =pr.get_eos_lnK(pvt,T,P,{'nitrogen':.83,'oxygen':.17},{'nitrogen':.75,'oxygen':.25})
##for i in lnK:
##    print i,'\t',math.exp(lnK[i])

##################### successive_substitution test#################
##print pr.successive_substitution(pvt,T,P,{'nitrogen':.83,'oxygen':.17},{'nitrogen':.75,'oxygen':.25})

##################### normalize comp####################
##print pr.normalise_comp({'nitrogen':83,'oxygen':127})

##################### stability incl tangent plane###############
##print pr.check_phase_stability(pvt,T,P)

##################### flash calculator ####################
##print pr.pvt_calculator(pvt,T,P)

##comp_l,comp_h,f= {'nitrogen':.83,'oxygen':.17},{'nitrogen':0.03,'oxygen':.97},.95
##comp_l,comp_h= pr.flash_itterate(pvt,T,P,comp_l,comp_h,f)

#################### re-arange matrix #####################
##A={1:{1:0,2:1,3:3,4:3},
##   2:{1:2,2:1,3:2,4:3},
##   3:{1:2,2:0,3:0,4:1},
##   4:{1:1,2:0,3:0,4:0},
##   }
##b={1:1,2:2,3:3,4:4}
##A,b= pr.re_arrange_square_matrix(A,b)
##A,b= pr.gauss_elim(A,b)
##
##for row in A:
##    string="|"
##    for col in A[row]:
##        string+=str(A[row][col])+"\t"
##    string+="|\t"+str(b[row])+" |"
##    print string

################################# Testing Matrix codes ############
##row,dimension=1,4
##row2=row
##if A[row][row]==0:
##        A_new,b_new={},{}
##        row3=row
##        while row3<=dimension:
##                b_new[row3-row+1]=b[row3]
##                col2=row
##                row_new={}
##                while col2<=dimension:
##                        row_new[col2-row+1]=A[row3][col2]
##                        col2+=1
##                A_new[row3-row+1]=row_new
##                row3+=1
##        A_new,b_new=pr.re_arrange_square_matrix(A_new,b_new)
##        row3=row
##        while row3<=dimension:
##                col2=row
##                b[row3]=b_new[row3-row+1]
##                while col2<=dimension:
##                        A[row3][col2]=A_new[row3-row+1][col2-row+1]
##                        col2+=1
##                row3+=1
##A,b=pr.gauss_elim(A,b)
##for row in A:
##    string="|"
##    for col in A[row]:
##        string+=str(A[row][col])+"\t"
##    string+="|\t"+str(b[row])+" |"
##    print string,'\n'
##print pr.solve_triangle_matrix(A,b)

############################## Isotherm plot ####################
##data=[]
##data2=[]
##data3=[]
##T0=int(-273.)
##T1=int(600-273)
##P0=int(1)
##P1=int(12000.6e3)
##steps=20
##for p in range(P0,P1,(P1-P0)/(steps)):
##    for T in range(T0,T1,(T1-T0)/steps):
##        V=None
##        V=pr.pvt_calculator(pvt,T,p)['V']
##
##        if 1e-4<V<.999:
##            data.append({'T':T,'P':p,'V':V})
##        elif V==1:
##            data2.append({'T':T,'P':p,'V':V})
##        elif V==0:
##            data3.append({'T':T,'P':p,'V':V})
##write.write_VLE(data)
##write.write_VLE2(data2)
##write.write_VLE3(data3)

##data=[]
##T0=int(-273.)
##T1=int(600-273)
##P0=int(1)
##P1=int(5e6)
##steps=100
##T=-50.
##for p in range(P0,P1,(P1-P0)/(steps)):
##    roots = pr.solve_PR_for_V(pvt,T,p)[0]
##    for current_root in range(0, len(roots)):
##        print(str(p) + "," + str(roots[current_root]))
        
##write.write_VLE(data)
##print pr.pvt_calculator(pvt,400-273,12e6)['V']
##print pr.pvt_calculator(pvt,400-273.15,90e5)
##
##print pr.solve_PR_for_V(pvt,300-273.15,5e5)[0]
##print pr.fug_minimum_gibbs(pvt,300-273.15,10.5e5,1)
##print pr.get_wilson_K(pvt,300-273.15,10.5e5)

##print pr.michelson_stability(pvt,T,P)
##print pr.check_phase_stability(pvt,T,P)
##print pr.get_initial_bubble_point(pvt)
# print pr.calculate_envelope(pvt)

# import flash
# import flash2
# print flash.get_binary_tangent_planes(pvt,T,P)
# import permutate
# print permutate.get_perms(2)


import pr_thermo
import TRAPP,chung_lee_starling,PVT
P=1.013e5
T=25
pvt= setPVT.get_props_and_equalize([(36,'methane'),(64,'water')])#,])

# print PVT.get_props(pvt,T,P)
# p=pr_thermo.Phase(pvt,{'type':'gas','frac':1.,'comp':{'water': 0.99966494548858955, 'methane': 1-0.99966494548858955}},T,P) # 0.9303327105473993
# p._set_density()
# print p.mol_density,p.mass_density,'\n'
# p._set_enthalpy()
# p._set_Cp()
# print p.mol_h
# print p.u_vap
# print TRAPP.get_visc_and_conductivity(p)
# print chung_lee_starling.get_visc_and_conductivity(p)
########################################### Test VLE #############################################################
# import flash3
# from matplotlib import pyplot as plt
# P=1.013e5
# T=25
# pvt= setPVT.get_props_and_equalize([(28,'npentane'),(15,'hexane'),(15,'heptane'),(14,'octane'),(14,'nonane'),(14,'decane')])
# x_1=[]
# y_1=[]
# c_1=[]
# x_2=[]
# y_2=[]
# c_2=[]
# x_3=[]
# y_3=[]
# c_3=[]
# for P in range(1,int(3226750.0),int(.5e5)):
#     for T in range(-273,600-273,10):
#         try:
#             a=flash3.Flash(pvt,T,P)
#             a.equilibriate()
#             if 0<a.phases[0]['frac']<1:
#
#                 # col='#%02x%02x%02x' % (a.phases[0]['frac']*255, 0, 255*(1-a.phases[0]['frac']))
#                 if a.phases[0]['frac']<.5:
#                     if a.phases[0]['frac']<.2:
#                         x_1.append(T+273)
#                         y_1.append(P / 1e3)
#                         col='#ff0000'
#                         c_1.append(col)
#                     else:
#                         x_2.append(T + 273)
#                         y_2.append(P / 1e3)
#                         col = '#00ff00'
#                         c_2.append(col)
#                 else:
#                     x_3.append(T + 273)
#                     y_3.append(P / 1e3)
#                     col='#0000ff'
#                     c_3.append(col)
#         except:
#             pass
# plt.scatter(x_1,y_1,c=c_1,alpha=1.,label='0 < V < 0.2')
# plt.scatter(x_2,y_2,c=c_2,alpha=1.,label='0.2 < V < 0.5')
# plt.scatter(x_3,y_3,c=c_3,alpha=1.,label='0.5 < V < 1')
# x=[340.877, 345.056, 350.056, 355.056, 360.056, 364.452, 365.056, 370.056, 375.056, 375.086, 380.056, 385.056, 386.526, 390.056, 395.056, 398.866, 400.056, 402.188, 405.056, 410.056, 412.208, 415.056, 420.056, 425.056, 426.663, 442.347, 459.374, 477.837, 492.644, 497.773, 519.055, 526.512, 540.988, 547.39, 551.193, 555.254, 555.259, 555.235, 555.081, 554.588, 554.216, 553.505, 552.748, 549.238, 543.493, 532.508, 525.87, 525.166, 516.9, 505.274, 492.63, 489.107, 467.032, 447.739, 443.269, 425.056, 420.056, 419.781, 415.056, 410.056, 405.056, 400.056, 398.791, 395.056, 390.056, 385.056, 380.056, 379.923, 375.056, 370.056, 365.056, 362.881, 360.056, 355.056, 350.056, 347.418, 345.056, 342.835, 340.056, 335.056, 333.33, 330.056, 325.056, 320.445, 320.056, 315.056, 310.056, 308.614, 305.056, 300.056, 297.715, 295.056, 290.056, 287.64, 285.056, 280.056, 278.297, 275.056, 270.056, 269.608, 0.0555556, 0.0555556, 0.0555556, 0.0555556, 0.0555556, 0.0555556, 0.0555556, 0.0555556, 0.0555556, 0.0555556, 0.0555556, 0.0555556, 0.0555556, 0.0555556, 0.0555556, 0.0555556, 0.0555556, 0.0555556, 0.0555556, 0.0555556, 0.0555556, 0.0555556, 0.0555556, 0.0555556, 0.0555556, 0.0555556, 0.0555556, 0.0555556, 0.0555556, 0.0555556, 0.0555556, 0.0555556, 0.0555556, 0.0555556, 0.0555556, 0.0555556, 0.0555556, 0.0555556, 0.0555556, 0.0555556, 0.0555556, 0.0555556, 0.0555556, 0.0555556, 0.0555556, 0.0555556, 0.0555556, 0.0555556, 0.0555556, 0.0555556]
# y=[10.0, 12.090399999999999, 15.0606, 18.6155, 22.840400000000002, 27.1829, 27.8287, 33.6815, 40.5078, 40.552099999999996, 48.4252, 57.559400000000004, 60.4966, 68.04480000000001, 80.0243, 90.25030000000001, 93.6497, 100.0, 109.082, 126.491, 134.638, 146.055, 167.965, 192.42, 200.856, 299.642, 447.013, 666.864, 900.002, 994.845, 1484.13, 1700.0, 2214.07, 2500.0, 2704.27, 3049.05, 3056.82, 3072.01, 3094.06, 3116.43, 3123.01, 3126.75, 3124.02, 3075.67, 2951.47, 2673.13, 2500.0, 2481.74, 2270.07, 1985.94, 1700.0, 1625.0, 1202.78, 900.002, 838.269, 617.24, 564.697, 561.909, 515.452, 469.389, 426.39, 386.338, 376.659, 349.116, 314.603, 282.68, 253.228, 252.482, 226.127, 201.26, 178.508, 169.244, 157.755, 138.887, 121.79, 113.447, 106.352, 100.0, 92.4655, 80.02289999999999, 76.04610000000001, 68.9205, 59.057199999999995, 50.975199999999994, 50.335300000000004, 42.660199999999996, 35.941199999999995, 34.1697, 30.0913, 25.0274, 22.9046, 20.6707, 16.9468, 15.353399999999999, 13.7857, 11.121799999999999, 10.2917, 8.89445, 7.04734, 6.89874, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
# plt.plot(x,y,label='PVTP')
# plt.legend(bbox_to_anchor=(0.05, .95), loc=2, borderaxespad=0.)
# plt.xlabel('T (K)')
# plt.ylabel('P (kPa)')
# plt.show()


#################################################################################################################
# print 'start'
# stuff=flash.System(pvt,P,T)
# vap=stuff._equilibriate()
# for i in stuff.phases:
#     print i['type'],i['phase'].moles,i['phase'].pvt['methane']['comp']


# print stuff.tangents
# print [(l['type'],l['phase'].moles) for l in stuff.phases]

# print pr.fug_minimum_gibbs(pvt,T, P, 1, phase='light',comp={'water':.0645,'propane':.9355})
# print pr.fug_minimum_gibbs(setPVT.get_props_and_equalize([(1e8,'water'),(1e-4,'aniline')]),T, P, 1, phase='light')