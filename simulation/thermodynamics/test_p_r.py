import sys
from os import path
sys.path.append(path.dirname( path.abspath(__file__)) )
import peng_robinson as pr
import setup_comp_data as setPVT
sys.path.append(path.dirname( path.dirname( path.dirname( path.abspath(__file__)))))
import data_writer as write
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

import flash

# print flash.get_binary_tangent_planes(pvt,T,P)
P=1.013e5
T=313.15-273.15
pvt= setPVT.get_props_and_equalize([(8,'water'),(2,'aniline')])
pvt['water']['k']={'aniline':0,'water':0}
pvt['aniline']['k']={'aniline':0,'water':0}
print pvt
stuff=flash.Prep(pvt,P,T)
print stuff.tangents
print stuff.pvt
print stuff.dominant
# print pr.fug_minimum_gibbs(pvt,T, P, 1, phase='light',comp={'water':.0645,'propane':.9355})
# print pr.fug_minimum_gibbs(setPVT.get_props_and_equalize([(1e8,'water'),(1e-4,'aniline')]),T, P, 1, phase='light')