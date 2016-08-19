from simulation.thermodynamics import setup_comp_data as set_comp
import flash3
import pr_thermo
from matplotlib import pyplot as plt
import copy


import peng_robinson


def get_props(pvt,T,P):
    # ################################################################################33
    # x_1=[]
    # y_1=[]
    # c_1=[]
    # z_1=[]
    # x_2=[]
    # y_2=[]
    # c_2=[]
    # x_3=[]
    # y_3=[]
    # c_3=[]
    # for P in range(1,int(60e5),int(1e5)):#range(1,int(45e5),int(1e5)):
    #     for T in range(160-270,800-273,5):#range(0,500,10):
    #         print 'something'
    #         try:
    #             print 'trying Flash'
    #             a=flash3.Flash(copy.deepcopy(pvt),T,P)
    #             print 'great success'
    #             a.equilibriate()
    #             if 0<=a.phases[0]['frac']<=1:
    #                 V=a.phases[0]['frac']
    #                 col='#%02x%02x%02x' % (V*255, 0, 255*(1-V))
    #                 if 0.01<V <.99 and P<3.2e6 and P>2:
    #                     x_1.append(T + 273)
    #                     y_1.append(P / 1e3)
    #                     z_1.append(V)
    #                     c_1.append(col)
                    # if a.phases[0]['frac']<.5:
                    #     if a.phases[0]['frac']<.2:
                    #         x_1.append(T+273)
                    #         y_1.append(P / 1e3)
                    #         col='#ff0000'
                    #         c_1.append(col)
                    #     else:
                    #         x_2.append(T + 273)
                    #         y_2.append(P / 1e3)
                    #         col = '#00ff00'
                    #         c_2.append(col)
                    # else:
                    #     x_3.append(T + 273)
                    #     y_3.append(P / 1e3)
                    #     col='#0000ff'
                    #     c_3.append(col)
            # except:
            #     pass
    # print '**************************************'
    # for i in range(len(x_1)):
    #     print x_1[i],y_1[i],z_1[i]
    # plt.scatter(x_1,y_1,c=c_1,alpha=1.,label='0 < V < 0.2',linewidths=0)
    # plt.scatter(x_2,y_2,c=c_2,alpha=1.,label='0.2 < V < 0.5')
    # plt.scatter(x_3,y_3,c=c_3,alpha=1.,label='0.5 < V < 1')
    # # x = ['321.5594444', '325.0455556', '330.0455556', '335.0455556', '340.0455556', '345.0455556', '345.0772222',
    # #      '350.0455556', '355.0455556', '355.6616667', '360.0455556', '365.0455556', '367.0272222', '370.0455556',
    # #      '375.0455556', '379.2566667', '380.0455556', '382.5433333', '385.0455556', '390.0455556', '392.4388889',
    # #      '395.0455556', '400.0455556', '405.0455556', '406.6655556', '410.0455556', '415.0455556', '420.0455556',
    # #      '422.0244444', '425.0455556', '438.5877778', '456.3883333', '470.5144444', '475.3711111', '495.2955556',
    # #      '502.18', '515.5066667', '521.5022222', '534.3133333', '534.3522222', '542.6644444', '547.2561111',
    # #      '547.3538889', '548.5205556', '548.4866667', '548.0233333', '547.05', '546.27', '545.0383333', '543.6261111',
    # #      '540.9088889', '540.3394444', '535.2172222', '530.1638889', '521.7683333', '514.7477778', '504.3466667',
    # #      '491.9366667', '482.8005556', '456.9866667', '452.385', '427.0927778', '425.0455556', '420.0455556', '415.995',
    # #      '415.0455556', '410.0455556', '405.0455556', '400.0455556', '395.0455556', '393.9688889', '390.0455556',
    # #      '385.0455556', '380.5133333', '380.0455556', '375.0455556', '370.0455556', '365.0455556', '360.0455556',
    # #      '359.3044444', '355.0455556', '350.0455556', '345.0455556', '342.7466667', '340.0455556', '335.0455556',
    # #      '330.0455556', '325.0455556', '323.4505556', '320.0455556', '319.4872222', '315.0455556', '310.0455556',
    # #      '305.0455556', '301.1486111', '300.0455556', '295.8359444', '295.0455556', '290.0455556', '285.4173889',
    # #      '285.0455556', '280.0455556', '275.0455556', '274.5075556', '271.2574444', '270.0455556', '265.0455556',
    # #      '260.0455556', '258.0504222', '257.6718222', '255.0455556', '250.0455556', '245.3469444', '245.0455556',
    # #      '243.9723333', '240.0455556', '235.0455556', '232.7426111', '232.4666667', '230.0455556', '225.0455556',
    # #      '222.5222222', '220.0455556', '219.7826111', '215.0455556', '213.7188889', '210.0455556', '205.8316667',
    # #      '205.7762222', '205.0455556', '200.0455556', '198.5055556', '195.0455556', '191.7766667', '190.0455556',
    # #      '189.7544444', '185.4988889', '185.0455556', '180.0455556', '179.6088889', '174.0577778', '168.8105556',
    # #      '168.6933333', '163.8394444', '159.1216667', '154.6383333', '150.3744444', '146.3155556', '142.4494444',
    # #      '138.765', '135.2511111', '131.8988889', '128.6988889', '125.6422222', '123.0455556', '122.7216667',
    # #      '119.9294444', '117.2577778', '114.7005556', '112.2511111', '109.9038889', '107.6533333', '105.4933333',
    # #      '103.4194444', '101.4272222', '99.51166667', '97.66888889', '95.895', '0.045555556', '0.045555556',
    # #      '0.045555556', '0.045555556', '0.045555556', '0.045555556', '0.045555556', '0.045555556', '0.045555556',
    # #      '0.045555556', '0.045555556', '0.045555556', '0.045555556', '0.045555556', '0.045555556', '0.045555556',
    # #      '0.045555556', '0.045555556', '0.045555556', '0.045555556', '0.045555556', '0.045555556', '0.045555556',
    # #      '0.045555556', '0.045555556', '0.045555556', '0.045555556', '0.045555556', '0.045555556', '0.045555556',
    # #      '0.045555556', '0.045555556', '0.045555556', '0.045555556', '0.045555556', '0.045555556', '0.045555556',
    # #      '0.045555556', '0.045555556', '0.045555556', '0.045555556', '0.045555556', '0.045555556', '0.045555556',
    # #      '0.045555556', '0.045555556', '0.045555556', '0.045555556', '0.045555556', '0.045555556']
    # # y = ['10.02773458', '11.75073436', '14.65273758', '18.13183196', '22.27489144', '27.17637419', '27.2101585',
    # #      '32.93977052', '39.67656863', '40.5796439', '47.50804732', '56.56500011', '60.52417643', '66.98787322',
    # #      '78.92724812', '90.27787973', '92.54432425', '100.0277556', '108.0110293', '125.5111912', '134.6652221',
    # #      '145.2400557', '167.4058031', '192.2306515', '200.8835715', '219.9496431', '250.8154017', '285.096823',
    # #      '299.6688919', '323.0786604', '447.0401859', '666.8919234', '900.0277893', '994.8720666', '1484.158497',
    # #      '1700.033339', '2214.092631', '2500.031993', '3300.030648', '3303.029867', '4100.036198', '4900.034852',
    # #      '4927.524248', '5588.283286', '5700.040402', '6018.488543', '6327.063393', '6500.039057', '6718.306379',
    # #      '6922.908293', '7241.563277', '7300.030816', '7747.224755', '8100.029471', '8572.734011', '8900.028126',
    # #      '9305.37089', '9700.026781', '9940.722747', '10437.42104', '10500.02544', '10711.69448', '10719.00292',
    # #      '10730.724', '10733.6198', '10733.48191', '10727.20768', '10711.90132', '10687.28704', '10653.29588',
    # #      '10644.74638', '10609.78997', '10556.70034', '10500.02544', '10493.75121', '10420.87362', '10337.86075',
    # #      '10244.57469', '10140.94649', '10124.67486', '10026.70037', '9901.836318', '9766.078553', '9700.026781',
    # #      '9619.427071', '9461.812926', '9293.029275', '9113.145065', '9053.367522', '8922.091348', '8900.028126',
    # #      '8719.937073', '8506.682239', '8282.395794', '8100.029471', '8047.28458', '7841.131346', '7801.555441',
    # #      '7545.691008', '7300.030816', '7279.967074', '7004.728374', '6720.616122', '6689.562137', '6500.039057',
    # #      '6428.554216', '6129.225236', '5823.5255', '5700.040402', '5676.522386', '5512.723644', '5198.074513',
    # #      '4900.034852', '4880.874323', '4812.616228', '4562.819182', '4245.474202', '4100.036198', '4082.654515',
    # #      '3930.645808', '3620.154216', '3465.677185', '3315.826536', '3300.030648', '3019.482986', '2942.413392',
    # #      '2732.861043', '2500.031993', '2497.025879', '2457.56029', '2195.042418', '2116.945505', '1946.555375',
    # #      '1792.174871', '1713.133377', '1700.033339', '1514.633323', '1495.56932', '1294.428574', '1277.646735',
    # #      '1075.595882', '903.6544314', '900.0277893', '757.6807047', '634.0418537', '529.5807694', '441.5381699',
    # #      '367.5105429', '305.4101561', '253.4278251', '210.005335', '173.8016555', '143.6718432', '118.6375321',
    # #      '100.0277556', '97.86896582', '80.66348583', '66.4280879', '54.66411561', '44.95291932', '36.94417643',
    # #      '30.34520451', '24.91282547', '20.44364398', '16.77011745', '13.75159284', '11.2729277', '9.23897438',
    # #      '0.027579028', '0.027579028', '0.027579028', '0.027579028', '0.027579028', '0.027579028', '0.027579028',
    # #      '0.027579028', '0.027579028', '0.027579028', '0.027579028', '0.027579028', '0.027579028', '0.027579028',
    # #      '0.027579028', '0.027579028', '0.027579028', '0.027579028', '0.027579028', '0.027579028', '0.027579028',
    # #      '0.027579028', '0.027579028', '0.027579028', '0.027579028', '0.027579028', '0.027579028', '0.027579028',
    # #      '0.027579028', '0.027579028', '0.027579028', '0.027579028', '0.027579028', '0.027579028', '0.027579028',
    # #      '0.027579028', '0.027579028', '0.027579028', '0.027579028', '0.027579028', '0.027579028', '0.027579028',
    # #      '0.027579028', '0.027579028', '0.027579028', '0.027579028', '0.027579028', '0.027579028', '0.027579028',
    # #      '0.027579028']
    # # plt.plot(x,y,label='PVTP')
    # plt.legend(bbox_to_anchor=(0.05, .95), loc=2, borderaxespad=0.)
    # plt.xlabel('T (K)')
    # plt.ylabel('P (kPa)')
    # plt.show()

    ###################################################################################
    # print [('%s: '%i,pvt[i]['comp']) for i in pvt]
    pvt1=dict(pvt.copy())
    sys=flash3.Flash(pvt1,T,P)
    phases=sys.equilibriate()
    sys=pr_thermo.System(pvt.copy(),phases,T,P)
    props={'enthalpy':0,'r_g':2,'m_g':.00002,'Cp_g':2000,'Cp_l':4000,'m_l':.002,'r_l':400}
    liq=[]
    c=0
    tot_mass=0
    m_l=0
    V_l=0
    Cp_l=0
    for p in sys.phases:
        if p.type=='gas':
            props['r_g']=p.mass_density
            props['m_g']=p.viscosity/1000
            props['Cp_g']=p.mass_Cp
        else:
            liq.append((p.mol_mass*p.mol_frac,p.viscosity/1000))
            m_l+=p.mol_mass*p.mol_frac
            V_l+=p.mol_frac/p.mol_density
            Cp_l+=p.mol_Cp*p.mol_frac
            print '!!!!!!!!!!!!! ',p.mass_density
        tot_mass+=p.mol_mass*p.mol_frac
        props['enthalpy']+=p.mol_frac*p.mol_h
        c+=1
    props['enthalpy']=props['enthalpy']/tot_mass
    props['Cp_l']=Cp_l/m_l if m_l>0 else 1000.
    props['r_l']=m_l/V_l if V_l>0 else 700.
    props['vapQ']=1-m_l/tot_mass
    for i in range(len(liq)):
        if i==0:
            viscosity=liq[i][1]
            mass=liq[i][0]
        else:
            viscosity=viscosity*(2*viscosity+liq[i][1]-2*(viscosity-liq[i][1])*(liq[i][0]/sum([liq[i][0],mass])))/\
                      (2*viscosity+liq[i][1]+(viscosity-liq[i][1])*(liq[i][0]/sum([liq[i][0],mass])))
            mass+=liq[i][0]
    props['m_l']=viscosity

    return props