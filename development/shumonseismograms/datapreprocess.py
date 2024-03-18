#Construct csv file to be read by MOOSE
##import packages
import numpy as np

##input parameter
xmin = -360
xmax = 360
ymin = 0
ymax = 0
numofnodes = 721
xcoords = np.linspace(xmin,xmax,numofnodes).flatten().tolist()
ycoords = [0]
numofevents = 18 #number of dynamic events 

##loop over events
for i in range(1,numofevents+1):

    #read files
    time_event_i = np.loadtxt("./events_data/time_event"+str(i)+".txt",delimiter=",")
    vx_event_i = np.loadtxt("./events_data/vx_event"+str(i)+".txt",delimiter=",")
    vy_event_i = np.loadtxt("./events_data/vy_event"+str(i)+".txt",delimiter=",")

    #flatten the data
    time = time_event_i.flatten().tolist()

    #write vx data file
    with open('./events_data_piecewisemultilinear/events_data_vx_num'+str(i)+'.txt','w') as f:

        f.write('AXIS X')
        f.write('\n')
        print("save AXIS X ...")
        for item in xcoords:
            # write each item on a new line
            f.write("%s " % item)
        f.write('\n')
        f.write('AXIS Y')
        f.write('\n')
        print("save AXIS Y ...")
        for item in ycoords:
            # write each item on a new line
            f.write("%s " % item)
        f.write('\n')
        f.write('AXIS T')
        f.write('\n')
        print("save AXIS T ...")
        for item in time:
            # write each item on a new line
            f.write("%s " % item)
        f.write('\n')
        f.write('DATA')
        f.write('\n')
        count = 0
        for index_i in range(np.shape(vx_event_i)[1]):
            vxdata = vx_event_i[:,index_i].flatten().tolist()
            count += 1
            print("save DATA ... ", round(count/np.shape(vx_event_i)[1]*100, 3), " %")
            for item in vxdata:
            # write each item on a new line
                f.write("%s\n" % item)
    
    #write vx data file
    with open('./events_data_piecewisemultilinear/events_data_vy_num'+str(i)+'.txt','w') as f:

        f.write('AXIS X')
        f.write('\n')
        print("save AXIS X ...")
        for item in xcoords:
            # write each item on a new line
            f.write("%s " % item)
        f.write('\n')
        f.write('AXIS Y')
        f.write('\n')
        print("save AXIS Y ...")
        for item in ycoords:
            # write each item on a new line
            f.write("%s " % item)
        f.write('\n')
        f.write('AXIS T')
        f.write('\n')
        print("save AXIS T ...")
        for item in time:
            # write each item on a new line
            f.write("%s " % item)
        f.write('\n')
        f.write('DATA')
        f.write('\n')
        count = 0
        for index_i in range(np.shape(vy_event_i)[1]):
            vydata = vy_event_i[:,index_i].flatten().tolist()
            count += 1
            print("save DATA ... ", round(count/np.shape(vy_event_i)[1]*100, 3), " %")
            for item in vydata:
            # write each item on a new line
                f.write("%s\n" % item)

    #dump time in csv file
    np.savetxt('./events_data_piecewisemultilinear/events_data_time_num'+str(i)+'.csv',time_event_i)

