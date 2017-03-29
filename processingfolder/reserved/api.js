$(function() {
    //activate home page and hide user information
    activaTab('homedesc');
    $('#user').hide();
    //Customize by setting base_url to cybercom/api docker application
    base_url = "/api";
    da_usertasks_url = base_url + "/queue/usertasks/.json?page_size=1000&taskname=ecopadq.tasks.tasks.teco_spruce_data_assimilation"
    //No other alterations is need to get the standard applicaiton running!
    login_url = base_url + "/api-auth/login/?next=";
    logout_url = base_url + "/api-auth/logout/?next=";
    user_task_url = base_url + "/queue/usertasks/.json?page_size=10";
    user_url = base_url + "/user/?format=json";
    prevlink=null;nextlink=null;forecast_task_id="default";
    
    $("#aprofile").click(function(){activaTab('profile')})
    $("#alogout").click(function(){window.location = logout_url.concat(document.URL);})
    //load_task_history(user_task_url);
    $('#prevlink').click(function(){load_task_history(prevlink);});
    $('#nextlink').click(function(){load_task_history(nextlink);});
    Handlebars.registerHelper('json_metatags', function(context) {
        if (typeof context !== 'undefined') {
            return JSON.stringify(context).replace(/"/g,'').replace(/\[/g,'').replace(/\]/g,'').replace(/,/g,', ');
        }else{
            return ""
        } 
    });
    Handlebars.registerHelper('time_zone',function(context){
        temp = new Date(context + "Z");
        return temp.toLocaleDateString() + "  " + temp.toLocaleTimeString();
    });
    //Default EcoPAD  Initial Parameters values 
    data= {"lat":47.50, "longitude":-93.46,"wsmax":60.1,"wsmin":0.2,"LAIMAX":5.3,"LAIMIN":0.3,"rdepth":150,"Rootmax":500,"Stemmax":1000,"SapR":1.0,"SapS":0.2,"SLA":40.0,"GLmax":39.2,"GRmax":20.25,"Gsmax":20.25,"stom_n":2,"a1":8,"Ds0":2000,"Vcmx0":80,"extkU":0.51,"xfang":0,"alpha":0.385,"Tau_Leaf":1.5,"Tau_Wood":40.0,"Tau_Root":0.8,"Tau_F":0.3,"Tau_C":5.86,"Tau_Micro":0.4,"Tau_SlowSOM":356.94,"Tau_Passive":2050.0,"gddonset":140.0,"Q10":2.0,"Rl0":30.2,"Rs0":7,"Rr0":29.0,"da_wsmax":0,"da_wsmin":0,"da_LAIMAX":0,"da_LAIMIN":0,"da_rdepth":0,"da_Rootmax":0,"da_Stemmax":0,"da_SapR":0,"da_SapS":0,"da_SLA":1,"da_GLmax":1,"da_GRmax":1,"da_Gsmax":1,"da_stom_n":0,"da_a1":0,"da_Ds0":0,"da_Vcmx0":1,"da_extkU":0,"da_xfang":0,"da_alpha":0,"da_Tau_Leaf":1,"da_Tau_Wood":1,"da_Tau_Root":1,"da_Tau_F":1,"da_Tau_C":1,"da_Tau_Micro":1,"da_Tau_SlowSOM":1,"da_Tau_Passive":1,"da_gddonset":1,"da_Q10":1,"da_Rl0":1,"da_Rs0":1,"da_Rr0":1,"min_wsmax":50.0,"min_wsmin":0.2,"min_LAIMAX":4.0,"min_LAIMIN":0.3,"min_rdepth":100,"min_Rootmax":400,"min_Stemmax":500,"min_SapR":1,"min_SapS":0.2,"min_SLA":10.0,"min_GLmax":10.0,"min_GRmax":10.,"min_Gsmax":10.,"min_stom_n":2,"min_a1":8,"min_Ds0":2000,"min_Vcmx0":14,"min_extkU":0.51,"min_xfang":0,"min_alpha":0.385,"min_Tau_Leaf":0.5,"min_Tau_Wood":5.0,"min_Tau_Root":0.3,"min_Tau_F":0.1,"min_Tau_C":1.0,"min_Tau_Micro":0.05,"min_Tau_SlowSOM":5.0,"min_Tau_Passive":500.0,"min_gddonset":100.0,"min_Q10":1.0,"min_Rl0":10.0,"min_Rs0":4.5,"min_Rr0":10.0,"max_wsmax":60.0,"max_wsmin":0.2,"max_LAIMAX":4.0,"max_LAIMIN":0.3,"max_rdepth":100,"max_Rootmax":500,"max_Stemmax":1000,"max_SapR":1,"max_SapS":0.2,"max_SLA":200.0,"max_GLmax":50.0,"max_GRmax":30.0,"max_Gsmax":30.0,"max_stom_n":2,"max_a1":8,"max_Ds0":2000,"max_Vcmx0":180,"max_extkU":0.51,"max_xfang":0,"max_alpha":0.385,"max_Tau_Leaf":3.0,"max_Tau_Wood":800.0,"max_Tau_Root":2.0,"max_Tau_F":0.5,"max_Tau_C":20.0,"max_Tau_Micro":0.5,"max_Tau_SlowSOM":1000.0,"max_Tau_Passive":4000.0,"max_gddonset":160.0,"max_Q10":4.0,"max_Rl0":45.0,"max_Rs0":10.5,"max_Rr0":45.0}
    $('#workflow_link').click(function(){ if (!($('#user').is(":visible"))) {setup_auth_workflow();} });
    $('#forecast_select').change(function(){show_images();});
    //$('#forecast_images').click(function(){ if (!($('#user').is(":visible"))) {show_images();} });
    load_frontpage();
    setTimeout( function(){ 
        // Turn gifs off at 1 minute 
        $(".eco-img").each(function(itm,e){
            $(e).attr('src',e.src.replace('.gif', '_stat.gif'));
        });
        }  , 60000 );
    $('#resetDA').click(function(){forecast_task_id="default";$('#visual_daf').val("Default Data Assimilation");});
   show_images() 
});//End of Document Ready

function load_frontpage(){
    //main description
    $('#homedesc').empty();
    workflow_template = Handlebars.templates['tmpl-ecopad-desc'];
    desc_data ={"forcing_uncertainty_label":"Forcing Uncertainty","forcing_uncertainty_url":'bower_components/img/forcing_uncertainty.gif',
                "add_data_label":"Added Data","add_data_url":'bower_components/img/add_data.gif'}
    $('#homedesc').append(workflow_template(desc_data));
    //images
    workflow_template = Handlebars.templates['tmpl-ecopad-plot-results'];
    $('#front_viz').empty();
    //er_forcast.png       foliage_forecast.png gpp_forecast.png     root_forecast.png    soil_forecast.png    wood_forecast.png
    img_data ={'zero_label':'GPP Forecast','zero_url':'bower_components/img/gpp_forecast.gif',
                'one_label':'ER Forecast','one_url':'bower_components/img/er_forecast.gif',
                'two_label':'Foliage Forecast','two_url':'bower_components/img/foliage_forecast.gif',
                'three_label':'Wood Forecast','three_url':'bower_components/img/wood_forecast.gif',
                'four_label':'Root Forecast','four_url':'bower_components/img/root_forecast.gif',
                'five_label':'Soil Forecast','five_url':'bower_components/img/soil_forecast.gif'}
    $('#front_viz').append(workflow_template(img_data));
    //toggle gif and stat
    $('#toggle').click(function(){
        var stp_animation=false;
        if ($(".eco-img")[0].src.endsWith("_stat.gif")){
            stp_animation= true;
        }
        $(".eco-img").each(function(itm,e){ 
            if ( stp_animation ) { 
                location.reload(); 
                //$(e).attr('src',e.src.replace('_stat.gif', '.gif')); 
            }else { 
                $(e).attr('src',e.src.replace('.gif', '_stat.gif'));
            } 
        });
    });
}
function show_images(){
	var tag_val= $('#forecast_select').val()
	console.log(tag_val)
	var url="/api/catalog/data/forecast/public/?query={'filter':{'tag':'"+tag_val+"'},'sort':[('timestamp',-1)]}&page_size=1"
	$.getJSON(url,function(data){
            er_url=data.results[0].result_url+"/plot/er_forecast_weekly.png";
            gpp_url=data.results[0].result_url+"/plot/gpp_forecast_weekly.png"
            //set local time
            ctimestamp=data.results[0].timestamp;
            temp = new Date(ctimestamp + "Z");
            exact_date=temp.toLocaleDateString() + "  " + temp.toLocaleTimeString()
	    templete_data={'er_url':er_url,'exact_date':exact_date,'gpp_url':gpp_url}	
	    $('#two_images').empty();
            workflow_template = Handlebars.templates['tmpl-forecast'];
            $('#two_images').append(workflow_template(templete_data));		
	})
}
function setup_auth_workflow(){
    set_auth(base_url,login_url);
    $('#user').show();
    load_task_history(user_task_url);
    showDA_tasks();
    $('.tr-result').hide() 
    //Load Inital Parameters
    load_workflow(data);
    $('#task').change(function(){task_change();});
    $('#setModelParameter').click(function(){showInitParameters();});
    $('#runModel').click(function(){ submitWorkflow();});
    $('.da_param').hide();
    $('.forecast').hide();
    //Error Checking for DA Parmeters
    $('.da_chkbx').click(function(event){
        da_num=0
        da_num = check_da_number()
        if ($(this).is(':checked')){
            console.log(da_num)
            if (da_num>=19){
                event.preventDefault();
                event.stopPropagation();
                alert("18 is the maximun number of DA Parameters. Please uncheck one parameter to allow for addition.")
                return false
            }else{
                $(this).attr('checked',true);
                $(this).val("1");
                return true;
            }
        }else{
            console.log(da_num)
            if (da_num==0){
                event.preventDefault();
                event.stopPropagation();
                alert("Must have at least one DA Parameter.")
                return false

            }else{
                $(this).attr('checked',false);
                $(this).val("0");
                return true;
            }
        } 
    });
    $('#selectDA').click(function(){$('#myDAModal').modal('show');})

}
function task_change(){
    if ($('#task').val()=="DA"){
        $('.da_param').show();
        $('.forecast').hide();
        $('#setModelParameter').html("Set DA Initial Parameters")
    }else if ( $('#task').val()=="forcasting") {
        $('.da_param').hide();
        $('.forecast').show();
        $('#setModelParameter').html("Set Initial Parameters")
    }else{
        $('.da_param').hide();
        $('.forecast').hide();
        $('#setModelParameter').html("Set Initial Parameters")
    }
}
function submitWorkflow(){
    //model_type
    tasks= {
            "simulation":"ecopadq.tasks.tasks.teco_spruce_simulation",
            "DA":"ecopadq.tasks.tasks.teco_spruce_data_assimilation",
            "forcasting":"ecopadq.tasks.tasks.teco_spruce_forecast",
            "FRC_sim":"ecopadq.tasks.tasks.frc_01_sim",
            "New_task_FRC":"ecopadq.tasks.tasks.new_task",
            "NCB_sim":"ecopadq.tasks.tasks.new_code_build_sim1",
            "NCB_sim_x":"ecopadq.tasks.tasks.new_code_build_x_sim"
                       }
    //mtype = $("#task").prop('selectedIndex'
    task_name = tasks[$("#task").val()]
    tags= ["TECO Spruce Simulation"]new_code_build_x
    params = $('#parameters').serializeObject()
    clean_params(params,$("#task").val());
    //setup task_data
    task_data = {"function": task_name,"queue": "celery","args":[params],"kwargs":{},"tags":[]};
    if (task_name=="ecopadq.tasks.tasks.teco_spruce_data_assimilation"){
        tags=["TECO Spruce Data Assimilation"]
        $.each($('.da_chkbx'),function(idx,ob){ if($(ob).is(':checkbox')){if ($(ob).is(':checked')){tags.push($(ob).attr('name').substring(3));}}})
    }
    if (task_name=="ecopadq.tasks.tasks.teco_spruce_forecast"){
        tags=["TECO Spruce Forecast"]
        tags.push("End_Date " + $('#forecastdate').val())
        nd = new Date($('#forecastdate').val());
        if(isNaN(nd)){
            alert("Please check Forecast End Date. Format YYYY-mm-dd")
            return;
        }
        fyear=nd.getFullYear();
        if (fyear>2024){
            alert("Please correct Forecast End Date. Maximum Date 2024-12-31")
            return;
        }
        console.log(fyear)
        start = new Date(fyear, 0, 0);
        diff = nd-start;
        oneDay = 1000 * 60 * 60 * 24;
        fday = Math.ceil(diff/oneDay);
        if (fday>=367){
            alert("Please correct Forecast End Date. Maximum Date 2024-12-31")
            return;
        }
        console.log(fday);
        //Warming Treatment
        temp_treatment=$('#forecasttemp').val();
        tags.push("Warming " + $('#forecasttemp').val())
        try{
            num=parseFloat(temp_treatment)
            if( num<0.0 || num>9.0){
                alert("Warming Parameter must be between 0-9 degree celsius")
                return;
            }
        }catch(e){
            alert("Warming Parameter must be a number!")
            return;
        }
        
        //Co2 Treatment
        co2_treatment =$('#forecastco2').val();
        tags.push("CO2 Adjustment " + co2_treatment)
        try{
            num=parseFloat(co2_treatment)
            if( num<380.0 || num>900.0){
                alert("CO2 Adjustment Parameter must be between 380-900 ppm")
                return;
            }
        }catch(e){
            alert("CO2 Adjustment Parameter must be a number!")
            return;
        }
        if (temp_treatment.indexOf('.')<0){
            temp_treatment = temp_treatment + '.0'
        }
        if (co2_treatment.indexOf('.')<0){
            co2_treatment = co2_treatment + '.0'
        }
        task_data.args = [params,fyear,fday ]
        //"da_task_id":forecast_task_id
        task_data.kwargs = {"temperature_treatment":temp_treatment,"co2_treatment":co2_treatment}
        if (forecast_task_id != "default"){
            task_data.kwargs["da_task_id"]=forecast_task_id;
        }
    }
    task_data.tags=tags
    url = "/api/queue/run/" + task_name + "/.json"

    $.postJSON(url,task_data ,function(data){
        work_flow_temp = Handlebars.templates['tmpl-workflow_reset'];
        $('#home').empty();
        $('#home').append(work_flow_temp({'model':'TECO Spruce '}));
        //$('#task_result').empty();
        //$('#task_result').append("<pre>" + JSON.stringify(data,null, 4) + "</pre>");
        $('#home').append('<iframe width="100%" frameborder="0" id="myIframe" src="history_result_meta.html?task_id=' + data.task_id + '" style="min-height:420px;"></iframe>');
        //$('#task_result').urlize();
        load_task_history(user_task_url);
    });

}
function isNumeric(n) {
  return !isNaN(parseFloat(n)) && isFinite(n);
}
function check_da_number(){
    da_num=0
    $.each($('.da_chkbx'),function(idx,ob){ if($(ob).is(':checkbox')){if ($(ob).is(':checked')){da_num +=1; }}})
    return da_num
}
function clean_params(params,task_type){
    $.each(data,function(im,value){
        if (!(im in params)){
            params[im]="0";
        }
    });
    if (!(task_type=="DA")){
        $.each(params,function(im,value){
            if (im.startsWith("da_") || im.startsWith("min_") || im.startsWith("max_")){
                delete params[im];
            };
        });
    };
};
function showDA_tasks(){
     $.getJSON(da_usertasks_url, function(data){
        $('#da_result_tbody').empty();
        tr_template = Handlebars.templates['tmpl-da-tr']
        $.each(data.results, function(i, item) {
            load_success_da_only(item,tr_template);
        });
    });
}
function load_success_da_only(item,tmpl){
    $.getJSON(item.result,function(data){
        if (data.result.status =="SUCCESS" && item.tags.length>0){
            temp=item.task_name.split('.')
            item['task_name']= temp[temp.length-1]
            item.timestamp = item.timestamp;
            $('#da_result_tbody').append(tmpl(item))
            $('.select-da-id').click(function(ob){setDAid($(this));});
        }
    });
}
function setDAid(ob){
    forecast_task_id= $(ob).val()
    $('#visual_daf').val("Task ID: " + $(ob).val())
}
$.postJSON = function(url, data, callback,fail) {
    return jQuery.ajax({
        'type': 'POST',
        'url': url,
        'contentType': 'application/json',
        'data': JSON.stringify(data),
        'dataType': 'json',
        'success': callback,
        'error':fail,
        'beforeSend':function(xhr, settings){
            xhr.setRequestHeader("X-CSRFToken", getCookie('csrftoken'));
        }
    });
}; 
function showInitParameters(){
    $("#myParameterModal").modal('show');
}
function load_workflow(data){
    $('#home').empty();
    workflow_template = Handlebars.templates['tmpl-ecopad'];
    $('#home').append(workflow_template({}));
    //Load Parameters into modal pop up
    $('#myParameterModalbody').empty();
    workflow_template = Handlebars.templates['tmpl-workflow'];
    $('#myParameterModalbody').append(workflow_template(data));
}
function submit_user(){
    console.log(user_url)
    $.post( user_url,$('#user_form').serializeObject(),function(data){
        data.csrftoken = getCookie('csrftoken')
        $('#profile').empty();
        //source = $('#user-template').html()
        //user_template = Handlebars.compile(source);
        user_template = Handlebars.templates['tmpl-user']
        $('#profile').append(user_template(data))
        $('#user_form').hide()
        $('#view_form').show()
        $('#reset_password').click(function(){$('#pass_form').toggle(!$('#pass_form').is(':visible'));});
    })
    .fail(function(){ alert("Error Occured on User Update.")});
    //$('#user_form').hide()
    //$('#view_form').show()
    //var formData = JSON.parse($("#user_form").serializeArray());
    //console.log(formData);
    return false;
}
function edit_user(){
    $('#user_form').show()
    $('#view_form').hide()
    return false;
}
function set_password(){
    pass = $('#pass_form').serializeObject()
    if (pass.password !== pass.password2){
        alert("Passwords were not identical")
        return false;

    }
    $.post( user_url,$('#pass_form').serializeObject(),function(data){
        $('#reset_password').click(function(){$('#pass_form').toggle(!$('#pass_form').is(':visible'));});
        alert(JSON.stringify(data))
    })
    .fail(function(){ alert("Error Occured on Password Reset.")});
    return false;
}
function set_auth(base_url,login_url){
    $.getJSON( base_url + "/user/.json",function(data){
        $('#user').html(data['username'].concat( ' <span class="caret"></span> '));
        $("#user").append($('<img style="border-radius:80px;">').attr("src",data['gravator_url'] + "?s=40&d=mm") );
        data.csrftoken = getCookie('csrftoken')
        //source = $('#user-template').html()
        //user_template = Handlebars.compile(source);
        $('#profile').empty()
        user_template = Handlebars.templates['tmpl-user']
        $('#profile').append(user_template(data))
        $('#user_form').hide()
        $('#view_form').show() 
        $('#reset_password').click(function(){$('#pass_form').toggle(!$('#pass_form').is(':visible'));});
    })
    .fail(function() {
        activaTab('homedesc')
        var slink = login_url.concat(document.URL);
        window.location = slink
    });
}
function activaTab(tab){
    $('a[href="#' + tab + '"]').tab('show')
};
function load_task_history(url){
    $.getJSON(url, function(data){
    prevlink = data.previous;
    nextlink = data.next;
    if (prevlink == null){$('#li_prevlink').addClass("disabled");} else {$('#li_prevlink').removeClass("disabled");};
    if (nextlink == null){$('#li_nextlink').addClass("disabled");} else {$('#li_nextlink').removeClass("disabled");};
    setTaskDisplay(data);
    //source = $('#tr-template').html();
    //tr_template = Handlebars.compile(source);
    tr_template = Handlebars.templates['tmpl-tr']
    $('#result_tbody').html("")//clear table
    $.each(data.results, function(i, item) {
        temp=item.task_name.split('.')
        item['task_name']= temp[temp.length-1]
        item.timestamp = item.timestamp; // .substring(0,19).replace('T',' ')
        $('#result_tbody').append(tr_template(item)) 
    });
    });
}
function setTaskDisplay(data){
    if (data.count <= data.meta.page_size){
        $('#task_count').text('Task 1 - ' + data.count +  ' Total ' + data.count );
    }else{
        rec_start = data.meta.page_size*data.meta.page - data.meta.page_size +1;
        rec_end = "";
        if(data.meta.page_size*data.meta.page >= data.count){
            rec_end = data.count;
        }else{
            rec_end = data.meta.page_size*data.meta.page;
        }   
        $('#task_count').text('Task ' + rec_start + ' - ' + rec_end  +  ' Total ' + data.count )
    }

}
function showdaResult(url,id){
    $.getJSON(url + ".json" , function(data){
        data_res = {"args":data.args,"result":data.result}
        json_data = JSON.stringify(data_res,null, 4);
        //$("#"+ id).html(json_data);
        //$("#" + id).urlize();
        //$("." + id).show();
        alert(json_data);
    });
}
function showResult(url){
    //myModalLabel -->title
    $.getJSON(url + ".json" , function(data){
        json_data = JSON.stringify(data,null, 4);
        $("#myModalbody").html(json_data);
        $("#myModalbody").urlize();
        $("#myModal").modal('show');
    });
}
jQuery.fn.urlize = function() {
    if (this.length > 0) {
        this.each(function(i, obj){
            // making links active
            var x = $(obj).html();
            var list = x.match( /\b(http:\/\/|www\.|http:\/\/www\.)[^ <]{2,200}\b/g );
            if (list) {
                for ( i = 0; i < list.length; i++ ) {
                    var prot = list[i].indexOf('http://') === 0 || list[i].indexOf('https://') === 0 ? '' : 'http://';
                    x = x.replace( list[i], "<a target='_blank' href='" + prot + list[i] + "'>"+ list[i] + "</a>" );
                }

            }
            $(obj).html(x);
        });
    }
};
$.fn.serializeObject = function()
{
    var o = {};
    var a = this.serializeArray();
    $.each(a, function() {
        if (o[this.name] !== undefined) {
            if (!o[this.name].push) {
                o[this.name] = [o[this.name]];
            }
            o[this.name].push(this.value || '');
        } else {
            o[this.name] = this.value || '';
        }
    });
    return o;
};
function getCookie(name) {
    var cookieValue = null;
    if (document.cookie && document.cookie != '') {
        var cookies = document.cookie.split(';');
        for (var i = 0; i < cookies.length; i++) {
            var cookie = jQuery.trim(cookies[i]);
            // Does this cookie string begin with the name we want?
            if (cookie.substring(0, name.length + 1) == (name + '=')) {
                cookieValue = decodeURIComponent(cookie.substring(name.length + 1));
                break;
            }
        }
    }
    return cookieValue;
}
