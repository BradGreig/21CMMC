import numpy
import math
from scipy import interpolate
from decimal import *
import string
import pickle
import os
import multiprocessing
from array import array

# These are values set in ANAL_PARAMS.H and HEAT_PARAMS.H. Note: These must be changed by hand if any values in the .H files are changed
Z_HEAT_MAX = 35.
ZPRIME_STEP_FACTOR = 1.02
DELTA_R_HII_FACTOR = 1.1

FRACT_FLOAT_ERR = 1e-3

TWOPLACES = Decimal(10) ** -2       # same as Decimal('0.01')
FOURPLACES = Decimal(10) ** -4       # same as Decimal('0.0001')
SIXPLACES = Decimal(10) ** -6       # same as Decimal('0.000001')

WalkerID1 = '1.000000'
WalkerID2 = '1.000000'

def worker(iterable,command_list):
    os.system(command_list[iterable])    


def CreateFilestrings_and_Commands(R_OPTION,redshifts,R_STEPS,last_filter_R,Tvir_STEPS,PL_INDEX_STEPS,IncludeAlpha,R_MFP_MIN,R_BUBBLE_MAX,
				Tvir_MIN,Tvir_MAX,EFF_FACTOR_PL_INDEX_MIN,EFF_FACTOR_PL_INDEX_MAX,ALPHA_PL):

	command_list = []	
	filenames = []	
	command_file_remove = []

	for i in range(len(redshifts)):

		if R_OPTION==1:
			for j in range(R_STEPS):		
				R_MFP = R_MFP_MIN + (R_BUBBLE_MAX - R_MFP_MIN)*float(j)/(float(R_STEPS) - 1.)

				for k in range(Tvir_STEPS):

					T_VIR = Tvir_MIN + (Tvir_MAX - Tvir_MIN)*float(k)/(float(Tvir_STEPS) - 1.)

					if IncludeAlpha is True:

						for l in range(PL_INDEX_STEPS):

							ALPHA_PL = EFF_FACTOR_PL_INDEX_MIN + (EFF_FACTOR_PL_INDEX_MAX - EFF_FACTOR_PL_INDEX_MIN)*float(l)/(float(PL_INDEX_STEPS) - 1.)

							command = "./Createfcoll_ionisation_LC %s %s %s %s %s %s 0"%(WalkerID1,WalkerID2,Decimal(repr(redshifts[i])).quantize(SIXPLACES),Decimal(repr(R_MFP)).quantize(SIXPLACES),Decimal(repr(T_VIR)).quantize(SIXPLACES),Decimal(repr(ALPHA_PL)).quantize(SIXPLACES))
							command_list.append(command)

							filenames.append("Box_fcoll_z%s_%s_%s_%s_%i.txt"%(Decimal(repr(redshifts[i])).quantize(SIXPLACES),Decimal(repr(R_MFP)).quantize(SIXPLACES),Decimal(repr(T_VIR)).quantize(SIXPLACES),Decimal(repr(ALPHA_PL)).quantize(SIXPLACES),0))

							command = "rm Box_fcoll_z%s_%s_%s_%s_%i.txt"%(Decimal(repr(redshifts[i])).quantize(SIXPLACES),Decimal(repr(R_MFP)).quantize(SIXPLACES),Decimal(repr(T_VIR)).quantize(SIXPLACES),Decimal(repr(ALPHA_PL)).quantize(SIXPLACES),0)
							command_file_remove.append(command)

					else:

						ALPHA_PL = Decimal(repr(float(ALPHA_PL))).quantize(SIXPLACES)

						command = "./Createfcoll_ionisation_LC %s %s %s %s %s %s 0"%(WalkerID1,WalkerID2,Decimal(repr(redshifts[i])).quantize(SIXPLACES),Decimal(repr(R_MFP)).quantize(SIXPLACES),Decimal(repr(T_VIR)).quantize(SIXPLACES),ALPHA_PL)
						print command
						command_list.append(command)

						filenames.append("Box_fcoll_z%s_%s_%s_%s_%i.txt"%(Decimal(repr(redshifts[i])).quantize(SIXPLACES),Decimal(repr(R_MFP)).quantize(SIXPLACES),Decimal(repr(T_VIR)).quantize(SIXPLACES),ALPHA_PL,0))

						command = "rm Box_fcoll_z%s_%s_%s_%s_%i.txt"%(Decimal(repr(redshifts[i])).quantize(SIXPLACES),Decimal(repr(R_MFP)).quantize(SIXPLACES),Decimal(repr(T_VIR)).quantize(SIXPLACES),ALPHA_PL,0)
						command_file_remove.append(command)

		if R_OPTION==0:

			R_MFP = last_filter_R
			for k in range(Tvir_STEPS):

				T_VIR = Tvir_MIN + (Tvir_MAX - Tvir_MIN)*float(k)/(float(Tvir_STEPS) - 1.)

				if IncludeAlpha is True:
			
					for l in range(PL_INDEX_STEPS):

						ALPHA_PL = EFF_FACTOR_PL_INDEX_MIN + (EFF_FACTOR_PL_INDEX_MAX - EFF_FACTOR_PL_INDEX_MIN)*float(l)/(float(PL_INDEX_STEPS) - 1.)

						command = "./Createfcoll_ionisation_LC %s %s %s %s %s %s 1"%(WalkerID1,WalkerID2,Decimal(repr(redshifts[i])).quantize(SIXPLACES),Decimal(repr(last_filter_R)).quantize(SIXPLACES),Decimal(repr(T_VIR)).quantize(SIXPLACES),Decimal(repr(ALPHA_PL)).quantize(SIXPLACES))
						command_list.append(command)

						filenames.append("Box_fcoll_z%s_%s_%s_%s_%i.txt"%(Decimal(repr(redshifts[i])).quantize(SIXPLACES),Decimal(repr(R_MFP)).quantize(SIXPLACES),Decimal(repr(T_VIR)).quantize(SIXPLACES),Decimal(repr(ALPHA_PL)).quantize(SIXPLACES),1))

						command = "rm Box_fcoll_z%s_%s_%s_%s_%i.txt"%(Decimal(repr(redshifts[i])).quantize(SIXPLACES),Decimal(repr(R_MFP)).quantize(SIXPLACES),Decimal(repr(T_VIR)).quantize(SIXPLACES),Decimal(repr(ALPHA_PL)).quantize(SIXPLACES),1)
						command_file_remove.append(command)

				else:

					ALPHA_PL = Decimal(repr(float(ALPHA_PL))).quantize(SIXPLACES)

					command = "./Createfcoll_ionisation_LC %s %s %s %s %s %s 1"%(WalkerID1,WalkerID2,Decimal(repr(redshifts[i])).quantize(SIXPLACES),Decimal(repr(last_filter_R)).quantize(SIXPLACES),Decimal(repr(T_VIR)).quantize(SIXPLACES),ALPHA_PL)
#					print command
					command_list.append(command)

					filenames.append("Box_fcoll_z%s_%s_%s_%s_%i.txt"%(Decimal(repr(redshifts[i])).quantize(SIXPLACES),Decimal(repr(R_MFP)).quantize(SIXPLACES),Decimal(repr(T_VIR)).quantize(SIXPLACES),ALPHA_PL,1))

					command = "rm Box_fcoll_z%s_%s_%s_%s_%i.txt"%(Decimal(repr(redshifts[i])).quantize(SIXPLACES),Decimal(repr(R_MFP)).quantize(SIXPLACES),Decimal(repr(T_VIR)).quantize(SIXPLACES),ALPHA_PL,1)
					command_file_remove.append(command)


	return command_list,filenames,command_file_remove

def CreateTable(num_processes,command_list,filenames,command_file_remove,R_OPTION,redshifts,R_STEPS,Tvir_STEPS,PL_INDEX_STEPS,IncludeAlpha,R_BUBBLE_MAX,last_filter_R,
				Tvir_MIN,Tvir_MAX,EFF_FACTOR_PL_INDEX_MIN,EFF_FACTOR_PL_INDEX_MAX,ALPHA_PL):

	if IncludeAlpha is True:

#		This is used for the filter step (R_BUBBLE_MAX). This will be interpolated over R_MFP, T_VIR and ALPHA_PL in the drive_21cm_streamlined.c
		f_coll_table_first_step = numpy.zeros(len(redshifts)*R_STEPS*Tvir_STEPS*PL_INDEX_STEPS)
#	    This is used for the final filter step (unfiltered density field). Will be fixed for all R_MFP. Therefore, this table is only interpolated over T_VIR and ALPHA_PL
		f_coll_table_last_step = numpy.zeros(len(redshifts)*Tvir_STEPS*PL_INDEX_STEPS)

		total_samples = len(redshifts)*Tvir_STEPS*PL_INDEX_STEPS*(1 + R_STEPS)

	else:

#		This is used for the filter step (R_BUBBLE_MAX). This will be interpolated over R_MFP and T_VIR in the drive_21cm_streamlined.c
		f_coll_table_first_step = numpy.zeros(len(redshifts)*R_STEPS*Tvir_STEPS)
#	    This is used for the final filter step (unfiltered density field). Will be fixed for all R_MFP. Therefore, this table is only interpolated over T_VIR
		f_coll_table_last_step = numpy.zeros(len(redshifts)*Tvir_STEPS)

		total_samples = len(redshifts)*Tvir_STEPS*(1 + R_STEPS)

	f_coll_vals = numpy.zeros(total_samples)

	num_divisions = int(numpy.floor(len(command_list)/num_processes))

	counter = 0

	for i in xrange(num_divisions):
		processes = []

		for ii in xrange(num_processes):
			p = multiprocessing.Process(target=worker, args=(ii + num_processes*i,command_list))
			p.start()
			processes.append(p)

		for p in processes:
			p.join()

		for ii in xrange(num_processes):
			f_coll_val = numpy.loadtxt('%s'%(filenames[counter]), usecols=(0,))			

			f_coll_vals[counter] = f_coll_val

			os.system(command_file_remove[counter])
			counter += 1

	remainder = len(command_list)%num_processes

	processes = []

	for ii in xrange(remainder):
		p = multiprocessing.Process(target=worker, args=(ii + num_divisions*num_processes,command_list))
		p.start()
		processes.append(p)

	for p in processes:
		p.join()

	for ii in xrange(remainder):
		f_coll_val = numpy.loadtxt('%s'%(filenames[counter]), usecols=(0,))			

		f_coll_vals[counter] = f_coll_val

		os.system(command_file_remove[counter])
		counter += 1





	counter = 0

	if R_OPTION == 1:
		for i in range(len(redshifts)):
			for j in range(R_STEPS):
				for k in range(Tvir_STEPS):

					if IncludeAlpha is True:
						for l in range(PL_INDEX_STEPS):

							f_coll_table_first_step[l + PL_INDEX_STEPS*( k + Tvir_STEPS*( j + R_STEPS*i ) )] = f_coll_vals[counter]

							counter += 1

					else:

						f_coll_table_first_step[k + Tvir_STEPS*( j + R_STEPS*i )] = f_coll_vals[counter]

						counter += 1

		if IncludeAlpha is True:
			output_file = open('Ionisation_fcoll_table_Rmax%s_Tmin%s_Tmax%s_PLmin%s_PLmax%s_%d_%sMpc'%(Decimal(repr(R_BUBBLE_MAX)).quantize(SIXPLACES),Decimal(repr(Tvir_MIN)).quantize(SIXPLACES),Decimal(repr(Tvir_MAX)).quantize(SIXPLACES),Decimal(repr(EFF_FACTOR_PL_INDEX_MIN)).quantize(SIXPLACES),Decimal(repr(EFF_FACTOR_PL_INDEX_MAX)).quantize(SIXPLACES),HII_DIM,BOX_LEN), 'wb')
		else:
			output_file = open('Ionisation_fcoll_table_Rmax%s_Tmin%s_Tmax%s_PL%s_%d_%sMpc'%(Decimal(repr(R_BUBBLE_MAX)).quantize(SIXPLACES),Decimal(repr(Tvir_MIN)).quantize(SIXPLACES),Decimal(repr(Tvir_MAX)).quantize(SIXPLACES),Decimal(repr(ALPHA_PL)).quantize(SIXPLACES),HII_DIM,BOX_LEN), 'wb')
		
		f_coll_table_first_step.tofile(output_file, format = 'float64')
		output_file.close()

	if R_OPTION == 0:

		for i in range(len(redshifts)):
			for k in range(Tvir_STEPS):

				if IncludeAlpha is True:
				
					for l in range(PL_INDEX_STEPS):

						f_coll_table_last_step[l + PL_INDEX_STEPS*( k + Tvir_STEPS*i )] = f_coll_vals[counter]

						counter += 1
				else:

					f_coll_table_last_step[k + Tvir_STEPS*i] = f_coll_vals[counter]

					counter += 1

		if IncludeAlpha is True:
			output_file = open('Ionisation_fcoll_table_final_Rmax%s_Tmin%s_Tmax%s_PLmin%s_PLmax%s_%d_%sMpc'%(Decimal(repr(R_BUBBLE_MAX)).quantize(SIXPLACES),Decimal(repr(Tvir_MIN)).quantize(SIXPLACES),Decimal(repr(Tvir_MAX)).quantize(SIXPLACES),Decimal(repr(EFF_FACTOR_PL_INDEX_MIN)).quantize(SIXPLACES),Decimal(repr(EFF_FACTOR_PL_INDEX_MAX)).quantize(SIXPLACES),HII_DIM,BOX_LEN), 'wb')
		else: 
			output_file = open('Ionisation_fcoll_table_final_Rmax%s_Tmin%s_Tmax%s_PL%s_%d_%sMpc'%(Decimal(repr(R_BUBBLE_MAX)).quantize(SIXPLACES),Decimal(repr(Tvir_MIN)).quantize(SIXPLACES),Decimal(repr(Tvir_MAX)).quantize(SIXPLACES),Decimal(repr(ALPHA_PL)).quantize(SIXPLACES),HII_DIM,BOX_LEN), 'wb')

		f_coll_table_last_step.tofile(output_file, format = 'float64')
		output_file.close()



if __name__ == '__main__':

	# If the full spin temperature computation is to be performed, a redshift must be provided to which to perform the evolution down to.
	TsCalc_z = 6.0

	redshifts = []

	z_prime = TsCalc_z*1.0001
        
	while (z_prime < Z_HEAT_MAX):
		z_prime = ((1.+z_prime)*ZPRIME_STEP_FACTOR - 1.)

	prev_z_prime = Z_HEAT_MAX
	z_prime = ((1.+z_prime)/ ZPRIME_STEP_FACTOR - 1.)
        
	while (z_prime > TsCalc_z):

		redshifts.append(z_prime)

		prev_z_prime = z_prime
		z_prime = ((1.+prev_z_prime) / ZPRIME_STEP_FACTOR - 1.)

	redshifts = numpy.array(redshifts)

#	L_FACTOR taken as set in ANAL_PARAMS.H (Note: These must be changed by hand if any values in the .H files are changed)
	L_FACTOR = 0.620350491

#	R_BUBBLE_MIN taken as set in ANAL_PARAMS.H
	R_BUBBLE_MIN = L_FACTOR
#	R_BUBBLE_MAX is the maximum allowed range as set in 21CMMC.py. Could make it a global value in 21CMMC.py to make it a referenceable value... 
	R_BUBBLE_MAX = 50.0
    
#	Minimum and maximum allowed range for Tvir as set in 21CMMC.py. Again, could make it a global value in 21CMMC.py to make it a referenceable value... 
	Tvir_MIN = 4.0
	Tvir_MAX = 6.0

	IncludeAlpha = False

	Fiducial_Alpha = 0.0

#	Minimum and maximum allowed range for Alpha as set in 21CMMC.py. Again, could make it a global value in 21CMMC.py to make it a referenceable value... 
	EFF_FACTOR_PL_INDEX_MIN = -2.0
	EFF_FACTOR_PL_INDEX_MAX = 2.0
    
	R_STEPS = 40
	Tvir_STEPS = 2
	PL_INDEX_STEPS = 0

	cell_length_factor = L_FACTOR

#	Box length as set in INIT_PARAMS.H (Note: These must be changed by hand if any values in the .H files are changed)
	BOX_LEN = 75
#	BOX_LEN = 300
#	BOX_LEN = 150

#	HII_DIM: number of voxels, set in INIT_PARAMS.H (Note: These must be changed by hand if any values in the .H files are changed)
	HII_DIM = 50
#	HII_DIM = 200
#	HII_DIM = 100
	
	R_MFP_MIN = max(R_BUBBLE_MIN, (cell_length_factor*BOX_LEN/float(HII_DIM)))
	last_filter_R = max(cell_length_factor*BOX_LEN/float(HII_DIM), R_BUBBLE_MIN)

	if IncludeAlpha is False:
		PL_INDEX_STEPS = 1


	num_processes = 8

	command_list = []	
	filenames = []	
	command_file_remove = []

	create_file = open("f_coll_lightcone_data_%s_%sMpc.txt"%(HII_DIM,BOX_LEN),"w")
	create_file.write("R_MFP_UB    %s\n"%(R_BUBBLE_MAX))
	create_file.write("X_RAY_TVIR_LB    %s\n"%(Tvir_MIN))
	create_file.write("X_RAY_TVIR_UB    %s\n"%(Tvir_MAX))
	create_file.write("ZETA_PL_LB    %s\n"%(EFF_FACTOR_PL_INDEX_MIN))
	create_file.write("ZETA_PL_UB    %s\n"%(EFF_FACTOR_PL_INDEX_MAX))
	create_file.write("R_MFP_STEPS    %s\n"%(R_STEPS))
	create_file.write("TVIR_STEPS    %s\n"%(Tvir_STEPS))
	create_file.write("PL_STEPS    %s\n"%(PL_INDEX_STEPS))	
	create_file.close()
	
	
	R_OPTION = 1

	command_list, filenames, command_file_remove = CreateFilestrings_and_Commands(R_OPTION,redshifts,R_STEPS,last_filter_R,Tvir_STEPS,PL_INDEX_STEPS,IncludeAlpha,R_MFP_MIN,R_BUBBLE_MAX,
				Tvir_MIN,Tvir_MAX,EFF_FACTOR_PL_INDEX_MIN,EFF_FACTOR_PL_INDEX_MAX,Fiducial_Alpha)


	CreateTable(num_processes,command_list,filenames,command_file_remove,R_OPTION,redshifts,R_STEPS,Tvir_STEPS,PL_INDEX_STEPS,IncludeAlpha,R_BUBBLE_MAX,last_filter_R,Tvir_MIN,Tvir_MAX,EFF_FACTOR_PL_INDEX_MIN,EFF_FACTOR_PL_INDEX_MAX,Fiducial_Alpha)
	
	R_OPTION = 0

	command_list, filenames, command_file_remove = CreateFilestrings_and_Commands(R_OPTION,redshifts,R_STEPS,last_filter_R,Tvir_STEPS,PL_INDEX_STEPS,IncludeAlpha,R_MFP_MIN,R_BUBBLE_MAX,
				Tvir_MIN,Tvir_MAX,EFF_FACTOR_PL_INDEX_MIN,EFF_FACTOR_PL_INDEX_MAX,Fiducial_Alpha)


	CreateTable(num_processes,command_list,filenames,command_file_remove,R_OPTION,redshifts,R_STEPS,Tvir_STEPS,PL_INDEX_STEPS,IncludeAlpha,R_BUBBLE_MAX,last_filter_R,Tvir_MIN,Tvir_MAX,EFF_FACTOR_PL_INDEX_MIN,EFF_FACTOR_PL_INDEX_MAX,Fiducial_Alpha)
	





