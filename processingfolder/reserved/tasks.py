from celery.task import task
from dockertask import docker_task
from subprocess import call,STDOUT
from jinja2 import Template
from shutil import copyfile, move
from glob import glob
import requests,os
from pymongo import MongoClient
from datetime import datetime
#Default base directory 
basedir="/data/static/"
spruce_data_folder="/data/local/spruce_data"
host= 'ecolab.cybercommons.org'
host_data_dir = os.environ["host_data_dir"] 
# "/home/ecopad/ecopad/data/static"

#print "hello-world"

#New Example task

@task()
def new_code_build_x_sim(pars): # ,model_type="0", da_params=None):
    """ Setup task convert parameters from html portal
	to file, and store the file in input folder.
	call teco_spruce_model.
    """
    task_id = str(new_code_build_x_sim.request.id)
    resultDir = setup_result_directory(task_id)
    #create param file 
    param_filename = create_template('SPRUCE_pars',pars,resultDir,check_params)
    #Run Spruce TECO code 
    host_data_resultDir = "{0}/static/ecopad_tasks/{1}".format(host_data_dir,task_id)
    host_data_dir_spruce_data="{0}/local/spruce_data".format(host_data_dir)	
    docker_opts = "-v {0}:/data:z -v {1}:/spruce_data:z".format(host_data_resultDir,host_data_dir_spruce_data)
    docker_cmd = "{0} {1} {2} {3} {4} {5}".format("/data/{0}".format(param_filename),"/spruce_data/SPRUCE_forcing.txt",
                                    "/spruce_data/SPRUCE_obs.txt",
                                    "/data", 0 , "/spruce_data/SPRUCE_da_pars.txt")
    result = docker_task(docker_name="new_code_build_x",docker_opts=docker_opts,docker_command=docker_cmd,id=task_id)
    #Run R Plots
    #os.makedirs("{0}/graphoutput".format(host_data_resultDir)) #make plot directory
    docker_opts = "-v {0}:/usr/local/src/myscripts/graphoutput:z ".format(host_data_resultDir)
    docker_cmd = None
    result = docker_task(docker_name="ecopad_r",docker_opts=docker_opts,docker_command=docker_cmd,id=task_id)
   

    #Clean up result Directory
#    clean_up(resultDir)
    #Create Report
    report_data ={'zero_label':'GPP','zero_url':'/ecopad_tasks/{0}/plot/{1}'.format(task_id,'gpp.png'),
                'one_label':'ER','one_url':'/ecopad_tasks/{0}/plot/{1}'.format(task_id,'er.png'),
                'two_label':'Foliage','two_url':'/ecopad_tasks/{0}/plot/{1}'.format(task_id,'foliage.png'),
                'three_label':'Wood','three_url':'/ecopad_tasks/{0}/plot/{1}'.format(task_id,'wood.png'),
                'four_label':'Root','four_url':'/ecopad_tasks/{0}/plot/{1}'.format(task_id,'root.png'),
                'five_label':'Soil','five_url':'/ecopad_tasks/{0}/plot/{1}'.format(task_id,'soil.png')}
    report_data['title']="SPRUCE Ecological Simulation Task Report"
    report_data['description']="Simulations of carbon fluxes and pool sizes for SPRUCE experiment based on user defined initial parameters."

    report = create_report('report',report_data,resultDir)
    result_url ="http://{0}/ecopad_tasks/{1}".format(result['host'],result['task_id'])
    #report_url = "http://{0}/ecopad_tasks/{1}/{2}".format(result['host'],result['task_id'],"report.htm")
    #{"report":report_url,"data":result_url}
    return result_url
  
#@task()
#def teco_spruce_simulation(pars): # ,model_type="0", da_params=None):
#    """ Setup task convert parameters from html portal
#	to file, and store the file in input folder.
#	call teco_spruce_model.
#    """
#    task_id = str(teco_spruce_simulation.request.id)
#    resultDir = setup_result_directory(task_id)
#    #create param file 
#    param_filename = create_template('SPRUCE_pars',pars,resultDir,check_params)
#    #Run Spruce TECO code 
#    host_data_resultDir = "{0}/static/ecopad_tasks/{1}".format(host_data_dir,task_id)
#    host_data_dir_spruce_data="{0}/local/spruce_data".format(host_data_dir)	
#    docker_opts = "-v {0}:/data:z -v {1}:/spruce_data:z".format(host_data_resultDir,host_data_dir_spruce_data)
#    docker_cmd = "{0} {1} {2} {3} {4} {5}".format("/data/{0}".format(param_filename),"/spruce_data/SPRUCE_forcing.txt",
#                                    "/spruce_data/SPRUCE_obs.txt",
#                                    "/data", 0 , "/spruce_data/SPRUCE_da_pars.txt")
#    result = docker_task(docker_name="teco_spruce",docker_opts=docker_opts,docker_command=docker_cmd,id=task_id)
#    #Run R Plots
#    #os.makedirs("{0}/graphoutput".format(host_data_resultDir)) #make plot directory
#    docker_opts = "-v {0}:/usr/local/src/myscripts/graphoutput:z ".format(host_data_resultDir)
#    docker_cmd = None
#    result = docker_task(docker_name="ecopad_r",docker_opts=docker_opts,docker_command=docker_cmd,id=task_id)
#  
#
#    #Clean up result Directory
#    clean_up(resultDir)
#    #Create Report
#    report_data ={'zero_label':'GPP','zero_url':'/ecopad_tasks/{0}/plot/{1}'.format(task_id,'gpp.png'),
#                'one_label':'ER','one_url':'/ecopad_tasks/{0}/plot/{1}'.format(task_id,'er.png'),
#                'two_label':'Foliage','two_url':'/ecopad_tasks/{0}/plot/{1}'.format(task_id,'foliage.png'),
#                'three_label':'Wood','three_url':'/ecopad_tasks/{0}/plot/{1}'.format(task_id,'wood.png'),
#                'four_label':'Root','four_url':'/ecopad_tasks/{0}/plot/{1}'.format(task_id,'root.png'),
#                'five_label':'Soil','five_url':'/ecopad_tasks/{0}/plot/{1}'.format(task_id,'soil.png')}
#    report_data['title']="SPRUCE Ecological Simulation Task Report"
#    report_data['description']="Simulations of carbon fluxes and pool sizes for SPRUCE experiment based on user defined initial parameters."
#
#    report = create_report('report',report_data,resultDir)
#    result_url ="http://{0}/ecopad_tasks/{1}".format(result['host'],result['task_id'])
#    #report_url = "http://{0}/ecopad_tasks/{1}/{2}".format(result['host'],result['task_id'],"report.htm")
#    #{"report":report_url,"data":result_url}
#    return result_url


@task()
def new_code_build_sim(pars): # ,model_type="0", da_params=None):
    """ Setup task convert parameters from html portal
	to file, and store the file in input folder.
	call teco_spruce_model.
    """
    task_id = str(new_code_build_sim.request.id)
    resultDir = setup_result_directory(task_id)
    #create param file 
    param_filename = create_template('SPRUCE_parameterfile',pars,resultDir,check_params)
    #Run Spruce TECO code 
    host_data_resultDir = "{0}/static/ecopad_tasks/{1}".format(host_data_dir,task_id)
    host_data_dir_spruce_data="{0}/local/example_code_mip_data".format(host_data_dir)	
    docker_opts = "-v {0}:/data:z -v {1}:/example_code_mip_data:z".format(host_data_resultDir,host_data_dir_spruce_data)
    docker_cmd = "{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11}".format("/data/{0}".format(param_filename),
        "example_code_mip_data/daily_soilt_2011-2014.txt",
        "example_code_mip/daily_watertable2011_2014.txt",
        "/example_code_mip_data/EM1forcing2011-2015.txt",
        "example_code_mip_data/obs_CH4_for_MEMCMC.txt",
        "/example_code_mip_data/obs_SPRUCE.txt",
        "example_code_mip_data/SPRUCE_hummock_toplayer.txt",
        "example_code_mip_data/SPRUCE_Snow_Depth_2011-2014.txt",
        "SPRUCE_soilt.txt","SPRUCE_Watr_Table_Level_2011-2014.txt",
        "/data",
         0)
    result = docker_task(docker_name="new_code_build",docker_opts=docker_opts,docker_command=docker_cmd,id=task_id)
    #Run R Plots
    #os.makedirs("{0}/graphoutput".format(host_data_resultDir)) #make plot directory
#    docker_opts = "-v {0}:/usr/local/src/myscripts/graphoutput:z ".format(host_data_resultDir)
#    docker_cmd = None
#    result = docker_task(docker_name="ecopad_r",docker_opts=docker_opts,docker_command=docker_cmd,id=task_id)
#   
#
#    #Clean up result Directory
#    clean_up(resultDir)
#    #Create Report
    report_data ={'zero_label':'GPP','zero_url':'/ecopad_tasks/{0}/plot/{1}'.format(task_id,'gpp.png'),
                'one_label':'ER','one_url':'/ecopad_tasks/{0}/plot/{1}'.format(task_id,'er.png'),
                'two_label':'Foliage','two_url':'/ecopad_tasks/{0}/plot/{1}'.format(task_id,'foliage.png'),
                'three_label':'Wood','three_url':'/ecopad_tasks/{0}/plot/{1}'.format(task_id,'wood.png'),
                'four_label':'Root','four_url':'/ecopad_tasks/{0}/plot/{1}'.format(task_id,'root.png'),
                'five_label':'Soil','five_url':'/ecopad_tasks/{0}/plot/{1}'.format(task_id,'soil.png')}
    report_data['title']="SPRUCE Ecological Simulation Task Report"
    report_data['description']="Simulations of carbon fluxes and pool sizes for SPRUCE experiment based on user defined initial parameters."
#
    report = create_report('report',report_data,resultDir)
    result_url ="http://{0}/ecopad_tasks/{1}".format(result['host'],result['task_id'])
    #report_url = "http://{0}/ecopad_tasks/{1}/{2}".format(result['host'],result['task_id'],"report.htm")
    #{"report":report_url,"data":result_url}
    return result_url

#
#
#W           W
# w    w    w
#  w  w w  w
#   w     w
#
#
#
#
#
#New_task from reserv
@task()
def new_task(pars): # ,model_type="0", da_params=None):
    """ Setup task convert parameters from html portal
	to file, and store the file in input folder.
	call teco_spruce_model.
    """
    host_data_dir = "/home/vova/ecopad/data"
    task_id = str(new_task.request.id)
    resultDir = setup_result_directory(task_id)
    #create param file 
    param_filename = create_template('f0_SPRUCE_parameterfile',pars,resultDir,check_params)
    #Run Spruce TECO code 
    host_data_resultDir = "{0}/static/ecopad_tasks/{1}".format(host_data_dir,task_id)
    host_data_dir_spruce_data="{0}/local/frc_01_input".format(host_data_dir)	
    docker_opts = "-v {0}:/data:z -v {1}:/frc_01_input:z".format(host_data_resultDir,host_data_dir_spruce_data)
    docker_cmd = "{0} {1}".format("/data/{0}".format(param_filename),
#COMMENT f0 Working
        "/frc_01_input/f1_EM1forcing2011-2015.txt") #,
#COMMENT f1 Working
      #  "/frc_01_input/f2_SPRUCE_Water_Table_Level_2011-2014.txt") #,
#COMMENT f2 is not Working ERROR: At line 4371 of file /source/frc_01.f90 (unit = 111, file = 'fort.111')
   #     "/frc_01_input/f3_SPRUCE_Snow_Depth_2011-2014.txt") #,
#COMMENT f3 Working
  #      "/frc_01_input/f4_obs_SPRUCE.txt",
#COMMENT f4 In process
   #    "/frc_01_input/f5_SPRUCE_hummock_toplayer.txt",
#COMMENT f5 In process
   #    "/frc_01_input/f6_SPRUCE_soilt.txt",
#COMMENT f6 In process
   #    "/frc_01_input/f7_daily_soilt_2011-2014.txt",
#COMMENT f7 In process
    #   "/frc_01_input/f8_daily_watertable_2011-2014.txt",
#COMMENT f8 In process
    #   "/frc_01_input/f9_obs_CH4_for_MEMCMC.txt")
#COMMENT f9 In process
#!!!!
#Need one more command to send output back (out of container before termination).
#!!!!
    result = docker_task(docker_name="frc_01",docker_opts=docker_opts,docker_command=docker_cmd,id=task_id)
    #Run R Plots
    #os.makedirs("{0}/graphoutput".format(host_data_resultDir)) #make plot directory
#    docker_opts = "-v {0}:/usr/local/src/myscripts/graphoutput:z ".format(host_data_resultDir)
#    docker_cmd = None
#    result = docker_task(docker_name="ecopad_r",docker_opts=docker_opts,docker_command=docker_cmd,id=task_id)
#   
#
#    #Clean up result Directory
#    clean_up(resultDir)
#    #Create Report
    report_data ={'zero_label':'GPP','zero_url':'/ecopad_tasks/{0}/plot/{1}'.format(task_id,'gpp.png'),
                'one_label':'ER','one_url':'/ecopad_tasks/{0}/plot/{1}'.format(task_id,'er.png'),
                'two_label':'Foliage','two_url':'/ecopad_tasks/{0}/plot/{1}'.format(task_id,'foliage.png'),
                'three_label':'Wood','three_url':'/ecopad_tasks/{0}/plot/{1}'.format(task_id,'wood.png'),
                'four_label':'Root','four_url':'/ecopad_tasks/{0}/plot/{1}'.format(task_id,'root.png'),
                'five_label':'Soil','five_url':'/ecopad_tasks/{0}/plot/{1}'.format(task_id,'soil.png')}
    report_data['title']="SPRUCE Ecological Simulation Task Report"
    report_data['description']="Simulations of carbon fluxes and pool sizes for SPRUCE experiment based on user defined initial parameters."
#
    report = create_report('report',report_data,resultDir)
    result_url ="http://{0}/ecopad_tasks/{1}".format(result['host'],result['task_id'])
    #report_url = "http://{0}/ecopad_tasks/{1}/{2}".format(result['host'],result['task_id'],"report.htm")
    #{"report":report_url,"data":result_url}
    return result_url




#Not brocken but archived (some problems with fortrun code, which there fixed)
#@task()
#def new_task(pars): # ,model_type="0", da_params=None):
#    """ Setup task convert parameters from html portal
#	to file, and store the file in input folder.
#	call teco_spruce_model.
#    """
#    task_id = str(new_task.request.id)
#    resultDir = setup_result_directory(task_id)
#    #create param file 
#    param_filename = create_template('f0_SPRUCE_parameterfile',pars,resultDir,check_params)
#    #Run Spruce TECO code 
#    host_data_resultDir = "{0}/static/ecopad_tasks/{1}".format(host_data_dir,task_id)
#    host_data_dir_spruce_data="{0}/local/frc_01_input".format(host_data_dir)	
#    docker_opts = "-v {0}:/data:z -v {1}:/frc_01_input:z".format(host_data_resultDir,host_data_dir_spruce_data)
#    docker_cmd = "{0}".format("/data/{0}".format(param_filename))#,
#COMMENT code "./input/{0}" works with call getarg(1,parafile) code.
#     #   "./input/f1_EM1forcing2011-2015.txt")#,

        #"/frc_01_input/f2_SPRUCE_Water_Table_Level_2011-2014.txt",

       # "/frc_01_input/f3_SPRUCE_Snow_Depth_2011-2014.txt",

       # "/frc_01_input/f4_obs_SPRUCE.txt",

        #"/frc_01_input/f5_SPRUCE_hummock_toplayer.txt",

       # "/frc_01_input/f6_SPRUCE_soilt.txt",

       # "/frc_01_input/f7_daily_soilt_2011-2014.txt",

        #"/frc_01_input/f8_daily_watertable_2011-2014.txt",

       # "/frc_01_input/f9_obs_CH4_for_MEMCMC.txt")

#
#End of My IN rpocess @data@
#



#
#Start of My checked docker_cmd
#
#COMMENT: this docker smd do not make any errors however output do not created for some reason.
#
#    docker_cmd = "{0} {1} {2} {3} {4} {5} {6} {7} {8} {9}".format("./input/{0}".format(param_filename),
#COMMENT code "./input/{0}" works with call getarg(1,parafile) code.
#        "./input/f1_EM1forcing2011-2015.txt",
#COMMENT code "./input/f1_EM1forcing2011-2015.txt") DO not give error with "call getarg(2,climatefile)"
#       "./input/f2_SPRUCE_Water_Table_Level_2011-2014.txt",
#COMMENT f2 same as f1
#       "./input/f3_SPRUCE_Snow_Depth_2011-2014.txt",
#COMMENT f3 same
#        "./input/f4_obs_SPRUCE.txt",
#COMMENT f4 same
#       "./input/f5_SPRUCE_hummock_toplayer.txt",
#COMMENT f5 same
#       "./input/f6_SPRUCE_soilt.txt",
#COMMENT f6 same
#       "./input/f7_daily_soilt_2011-2014.txt",
#COMMENT f7 same
#       "./input/f8_daily_watertable_2011-2014.txt",
#COMMENT f8 same
#       "./input/f9_obs_CH4_for_MEMCMC.txt")
#COMMENT f9 same
#
#End of My checked docker_cmd
#

#
#Start of My original docker_cmd:
#
#    docker_cmd = "{0} {1} {2} {3} {4} {5} {6} {7} {8} {9}".format("/data/{0}".format(param_filename),
#        "/frc_01_input/f1_EM1forcing2011-2015.txt",
#        "/frc_01_input/f2_SPRUCE_Water_Table_Level_2011-2014.txt",
#        "/frc_01_input/f3_SPRUCE_Snow_Depth_2011-2014",
#        "/frc_01_input/f4_obs_SPRUCE.txt",
#        "/frc_01_input/f5_SPRUCE_hummock_toplayer.txt",
#        "/frc_01_input/f6_SPRUCE_soilt.txt",
#        "/frc_01_input/f7_daily_soilt_2011-2014.txt",
#        "/frc_01_input/f8_daily_watertable_2011-2014.txt",
#        "/frc_01_input/f9_obs_CH4_for_MEMCMC.txt")
#
#End of My original docker_cmd
#
#    result = docker_task(docker_name="frc_01",docker_opts=docker_opts,docker_command=docker_cmd,id=task_id)
#    #Run R Plots
#    #os.makedirs("{0}/graphoutput".format(host_data_resultDir)) #make plot directory
#    docker_opts = "-v {0}:/usr/local/src/myscripts/graphoutput:z ".format(host_data_resultDir)
#    docker_cmd = None
#    result = docker_task(docker_name="ecopad_r",docker_opts=docker_opts,docker_command=docker_cmd,id=task_id)
#   
#
#    #Clean up result Directory
#    clean_up(resultDir)
#    #Create Report
#    report_data ={'zero_label':'GPP','zero_url':'/ecopad_tasks/{0}/plot/{1}'.format(task_id,'gpp.png'),
#                'one_label':'ER','one_url':'/ecopad_tasks/{0}/plot/{1}'.format(task_id,'er.png'),
#                'two_label':'Foliage','two_url':'/ecopad_tasks/{0}/plot/{1}'.format(task_id,'foliage.png'),
#                'three_label':'Wood','three_url':'/ecopad_tasks/{0}/plot/{1}'.format(task_id,'wood.png'),
#                'four_label':'Root','four_url':'/ecopad_tasks/{0}/plot/{1}'.format(task_id,'root.png'),
#                'five_label':'Soil','five_url':'/ecopad_tasks/{0}/plot/{1}'.format(task_id,'soil.png')}
#    report_data['title']="SPRUCE Ecological Simulation Task Report"
#    report_data['description']="Simulations of carbon fluxes and pool sizes for SPRUCE experiment based on user defined initial parameters."
#
#    report = create_report('report',report_data,resultDir)
#    result_url ="http://{0}/ecopad_tasks/{1}".format(result['host'],result['task_id'])
#    #report_url = "http://{0}/ecopad_tasks/{1}/{2}".format(result['host'],result['task_id'],"report.htm")
#    #{"report":report_url,"data":result_url}
#    return result_url

def clean_up(resultDir):

    move("{0}/SPRUCE_pars.txt".format(resultDir),"{0}/input/SPRUCE_pars.txt".format(resultDir))
    move("{0}/SPRUCE_yearly.txt".format(resultDir),"{0}/output/SPRUCE_yearly.txt".format(resultDir))
    for mvfile in glob("{0}/Simu_dailyflux*.txt".format(resultDir)):
        move(mvfile, "{0}/output".format(resultDir))
    for mvfile in glob("{0}/*.png".format(resultDir)):
        move(mvfile, "{0}/plot".format(resultDir))

    # Yuanyuan add to clear up forecast_csv
    #current_date=datetime.now().strftime("%Y-%m-%d")
    current_date=datetime.now().strftime("%Y")
    #if not os.path.exists("{0}/forecast_csv/{1}".format(basedir,current_date)):
    #    os.makedirs("{0}/forecast_csv/{1}".format(basedir,current_date))
    
    # make one folder for all the time, changed 01_04_2017
    if not os.path.exists("{0}/forecast_csv/ecopad_vdv".format(basedir)):
        os.makedirs("{0}/forecast_csv/ecopad_vdv".format(basedir))

    #for afile in glob.iglob("{0}/forecast_csv/{1}*".format(basedir,current_date)):
    #	print afile
    # 	os.remove(afile)
   
    try: 
        for mvfile in glob("{0}/*.csv".format(resultDir)):
            head,tail=os.path.split(mvfile)
            #dst_file=os.path.join("{0}/forecast_csv/{1}/{2}".format(basedir,current_date,tail))
            # modified 01_04_2017
            dst_file=os.path.join("{0}/forecast_csv/ecopad_vdv/{1}".format(basedir,tail))
            i=1 
            if os.path.exists(dst_file):
                with open(dst_file, 'a') as singleFile:
                    for line in open(mvfile, 'r'):
                       if i > 1:
                          singleFile.write(line)          
                          #print i
                       i=2
                os.remove(mvfile)
        else: 
            #move(mvfile,"{0}/forecast_csv/{1}".format(basedir,current_date))
            move(mvfile,"{0}/forecast_csv/ecopad_vdv".format(basedir)) 
    except:
        pass 

    try:
        move("{0}/SPRUCE_da_pars.txt".format(resultDir),"{0}/input/SPRUCE_da_pars.txt".format(resultDir))
        move("{0}/Paraest.txt".format(resultDir),"{0}/input/Paraest.txt".format(resultDir))
    except:
        pass

def create_template(tmpl_name,params,resultDir,check_function):
    tmpl = os.path.join(os.path.dirname(__file__),'templates/{0}.tmpl'.format(tmpl_name))
    with open(tmpl,'r') as f:
        template=Template(f.read())
    params_file = os.path.join(resultDir,'{0}.txt'.format(tmpl_name))
    with open(params_file,'w') as f2:
        f2.write(template.render(check_function(params)))
    return '{0}.txt'.format(tmpl_name)

def create_report(tmpl_name,data,resultDir):
    tmpl = os.path.join(os.path.dirname(__file__),'templates/{0}.tmpl'.format(tmpl_name))
    with open(tmpl,'r') as f:
        template=Template(f.read())
    report_file = os.path.join(resultDir,'{0}.htm'.format(tmpl_name))
    with open(report_file,'w') as f2:
        f2.write(template.render(data))
    return '{0}.htm'.format(tmpl_name)

def setup_result_directory(task_id):
    resultDir = os.path.join(basedir, 'ecopad_tasks/', task_id)
    os.makedirs(resultDir)
    os.makedirs("{0}/input".format(resultDir))
    os.makedirs("{0}/output".format(resultDir))
    os.makedirs("{0}/plot".format(resultDir))
    return resultDir 

def check_params(pars):
    """ Check params and make floats."""
    for param in ["latitude","longitude","wsmax","wsmin","LAIMAX","LAIMIN","SapS","SLA","GLmax","GRmax","Gsmax",
                    "extkU","alpha","Tau_Leaf","Tau_Wood","Tau_Root","Tau_F","Tau_C","Tau_Micro","Tau_SlowSOM",
                    "gddonset","Rl0" ]:
        try:
            inside_check(pars,param)
        except:
            pass
        try:
            inside_check(pars, "min_{0}".format(param))
        except:
            pass
        try:
            inside_check(pars, "max_{0}".format(param))
        except:
            pass
    return pars  

def inside_check(pars,param):
    if not "." in str(pars[param]):
        pars[param]="%s." % (str(pars[param]))
    else:
        pars[param]=str(pars[param])  
