

server = function(input, output){
  # pixel height for each panel. i.e row height when plotting by species
  height_per_panel = 350

  ref_table_reduced = summary_dt %>%
                as.data.frame(.)

  output$summarytable = DT::renderDataTable({
    summary_df = summary_dt %>%
                 as.data.frame(.,stringsAsFactors=FALSE)
    summary_DT = DT::datatable(summary_df, filter = 'top',rownames=FALSE,
    options = list(scrollX = TRUE, search = list(regex = TRUE, caseInsensitive = FALSE),pageLength = 25))
    return(summary_DT)
  })
  outputOptions(output, "summarytable", suspendWhenHidden = FALSE)

#   filtered_table = reactive({
#     req(input$summarytable_rows_selected)
#     keep_models = c(ref_table_reduced[input$summarytable_rows_selected, ]$id)
#     return(as.data.frame(summary_dt[id%in%keep_models],stringsAsFactors=FALSE))  
#   })

  filtered_id = reactive({
    req(input$summarytable_rows_selected)
    keep_models = c(ref_table_reduced[input$summarytable_rows_selected, ]$id)
    return(keep_models)  
  })

  # define plots
  output$ts_plots = renderPlot({
    input_models = unique(filtered_id())
    input_category = input$category
    if(length(input_models) < 1 | length(input_category) < 1){
      return()
    }
    
    selected_models = sapply(input_models,function(x)all_dirs[grep(x,all_dirs,fixed=TRUE)])
    plot_dt.list = list("SSB"=NA, "F"=NA, "Recruitment"=NA,"Depletion (D)"=NA,"SSB (Unfished)"=NA,"SSB/SSBmsy"=NA,"F/Fmsy"=NA)
    
    if("SSB" %in% input_category){
        plot_dt.list[["SSB"]] = rbindlist(lapply(selected_models,function(x)fread(file.path(model_stem,x,"ssb.csv")))) %>%
             setnames(.,c("Value","StdDev"),c("value","sd")) %>% 
             .[,cv:=sd/value] %>%
			 .[,l95:=exp(log(value)-2*sqrt(log(cv^2+1)))] %>%
			 .[,u95:=exp(log(value)+2*sqrt(log(cv^2+1)))] %>%
             .[,.(id,type,yr,value,l95,u95)] %>%
             .[,type:="SSB"] %>%
             .[,type:=factor(type,levels=c("SSB", "F", "Recruitment","Depletion (D)","SSB (Unfished)","SSB/SSBmsy","F/Fmsy"))]
    }

    if("F" %in% input_category){
        plot_dt.list[["F"]] = rbindlist(lapply(selected_models,function(x)fread(file.path(model_stem,x,"f.csv")))) %>%
             setnames(.,c("Value","StdDev"),c("value","sd")) %>% 
             .[,cv:=sd/value] %>%
			 .[,l95:=exp(log(value)-2*sqrt(log(cv^2+1)))] %>%
			 .[,u95:=exp(log(value)+2*sqrt(log(cv^2+1)))] %>%
             .[,.(id,type,yr,value,l95,u95)] %>%
             .[,type:="F"] %>%
             .[,type:=factor(type,levels=c("SSB", "F", "Recruitment","Depletion (D)","SSB (Unfished)","SSB/SSBmsy","F/Fmsy"))]
    }

    if("Recruitment" %in% input_category){
        plot_dt.list[["Recruitment"]] = rbindlist(lapply(selected_models,function(x)fread(file.path(model_stem,x,"rec.csv")))) %>%
             setnames(.,c("Value","StdDev"),c("value","sd")) %>% 
             .[,cv:=sd/value] %>%
			 .[,l95:=exp(log(value)-2*sqrt(log(cv^2+1)))] %>%
			 .[,u95:=exp(log(value)+2*sqrt(log(cv^2+1)))] %>%
             .[,.(id,type,yr,value,l95,u95)] %>%
             .[,type:="Recruitment"] %>%
             .[,type:=factor(type,levels=c("SSB", "F", "Recruitment","Depletion (D)","SSB (Unfished)","SSB/SSBmsy","F/Fmsy"))]
    }

    if("Depletion (D)" %in% input_category){
        plot_dt.list[["Depletion (D)"]] = rbindlist(lapply(selected_models,function(x)fread(file.path(model_stem,x,"dyn_dep.csv")))) %>%
             .[Era %in% c("TIME")] %>%
             setnames(.,c("Yr"),c("yr")) %>% 
             .[,value:=SSB/SSB_nofishing] %>%
			 .[,l95:=value] %>%
			 .[,u95:=value] %>%
             .[,type:="Depletion (D)"] %>%
             .[,.(id,type,yr,value,l95,u95)] %>%
             .[,type:=factor(type,levels=c("SSB", "F", "Recruitment","Depletion (D)","SSB (Unfished)","SSB/SSBmsy","F/Fmsy"))]
    }

    if("SSB (Unfished)" %in% input_category){
        plot_dt.list[["SSB (Unfished)"]] = rbindlist(lapply(selected_models,function(x)fread(file.path(model_stem,x,"dyn_dep.csv")))) %>%
             .[Era %in% c("TIME")] %>%
             setnames(.,c("Yr"),c("yr")) %>% 
             .[,value:=SSB_nofishing] %>%
			 .[,l95:=value] %>%
			 .[,u95:=value] %>%
             .[,type:="SSB (Unfished)"] %>%
             .[,.(id,type,yr,value,l95,u95)] %>%
             .[,type:=factor(type,levels=c("SSB", "F", "Recruitment","Depletion (D)","SSB (Unfished)","SSB/SSBmsy","F/Fmsy"))]
    }
    if("SSB/SSBmsy" %in% input_category){
        plot_dt.list[["SSB/SSBmsy"]] = rbindlist(lapply(selected_models,function(x)fread(file.path(model_stem,x,"kobe.csv")))) %>%
             setnames(.,c("Year"),c("yr")) %>% 
             .[,value:=B.Bmsy] %>%
			 .[,l95:=value] %>%
			 .[,u95:=value] %>%
             .[,type:="SSB/SSBmsy"] %>%
             .[,.(id,type,yr,value,l95,u95)] %>%
             .[,type:=factor(type,levels=c("SSB", "F", "Recruitment","Depletion (D)","SSB (Unfished)","SSB/SSBmsy","F/Fmsy"))]
    }
    if("F/Fmsy" %in% input_category){
        plot_dt.list[["F/Fmsy"]] = rbindlist(lapply(selected_models,function(x)fread(file.path(model_stem,x,"kobe.csv")))) %>%
             setnames(.,c("Year"),c("yr")) %>% 
             .[,value:=F.Fmsy] %>%
			 .[,l95:=value] %>%
			 .[,u95:=value] %>%
             .[,type:="F/Fmsy"] %>%
             .[,.(id,type,yr,value,l95,u95)] %>%
             .[,type:=factor(type,levels=c("SSB", "F", "Recruitment","Depletion (D)","SSB (Unfished)","SSB/SSBmsy","F/Fmsy"))]
    }

    plot_dt = rbindlist(plot_dt.list[which(sapply(plot_dt.list,length)>1)]) %>% .[order(id,type,yr)]

    if(nrow(plot_dt) == 0){
      return()
    }
    if(length(input_category)>2)
    {
      n_col = 3
    } else {
      n_col = length(input_category)
    }
    p = plot_dt %>%
      ggplot() +
			ylim(0,NA) +
			xlab("Year") +
			facet_wrap(~type,scales="free_y",ncol=n_col)

      p = p + ylab("Metric") +
              geom_hline(yintercept=0)
      if(input$se)
      {
        p = p + geom_ribbon(aes(x=yr,ymin=l95,ymax=u95,group=id,fill=id),alpha=0.25)
      }
	  p = p + geom_path(aes(x=yr,y=value,group=id,color=id),linewidth=1.5)
    
      p = p + viridis::scale_color_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			viridis::scale_fill_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
							panel.grid.major = element_line(color = 'gray70',linetype = "dotted"), 
							panel.grid.minor = element_line(color = 'gray70',linetype = "dotted"),
							strip.background =element_rect(fill="white"),
							legend.key = element_rect(fill = "white"))
			
    return(p)
  },
  height=function(){
    if(length(input$category)>2)
    {
      n_col = 3
    } else {
      n_col = length(input$category)
    }
    return(max(height_per_panel*1.5, (height_per_panel * ceiling(length(input$category) / n_col))))
  })
  # growth_plot
  output$growth_plot = renderPlot({
    input_models = unique(filtered_id())
    input_growth_var = input$growth_var
    input_bio_sex = input$bio_sex
    if(length(input_models) < 1|length(input_bio_sex)<1){
      return()
    }
    selected_models = sapply(input_models,function(x)all_dirs[grep(x,all_dirs,fixed=TRUE)])
    
    plot_dt = rbindlist(lapply(selected_models,function(x)fread(file.path(model_stem,x,"bio.csv")))) %>%
             .[,.(id,Sex,A1,A2,L_a_A1,L_a_A2,K,CVmin,CVmax)]
    max_age_dt = rbindlist(lapply(selected_models,function(x)fread(file.path(model_stem,x,"selex_a.csv"))[,.(max_age=max(as.numeric(variable))),by=id]))
    ages_dt.list = as.list(rep(NA,length(selected_models)))
    for(i in 1:length(selected_models)){
        ages_dt.list[[i]] = as.data.table(expand.grid(id=input_models[i],age=0:max_age_dt[id==input_models[i]]$max_age))
    }
    age_dt = rbindlist(ages_dt.list)

    plot_dt = merge(plot_dt,age_dt,by="id",allow.cartesian=TRUE) %>%
              .[,pred_len := L_a_A1 + (L_a_A2-L_a_A1)*((1-exp(-K*(age-A1)))/(1-exp(-K*(A2-A1))))] %>%
              .[age==A1,sd:=pred_len*CVmin] %>%
              .[age<A2,sd:=pred_len*(CVmin+((pred_len-L_a_A1)/(L_a_A2-L_a_A1))*(CVmax-CVmin))] %>%
              .[age==A2,sd:=pred_len*CVmax] %>%
              .[,u95:=pred_len+1.96*sd] %>%
              .[,l95:=pred_len-1.96*sd] %>%
              .[,.(id,Sex,age,pred_len,u95,l95)] %>%
              .[,Sex:=factor(Sex,levels=as.character(2:1),labels=c("Male","Female"))] %>%
              .[Sex %in% input_bio_sex]

    if(nrow(plot_dt) == 0){
      return()
    }
    
    n_col = length(input_bio_sex)
    
    p = plot_dt %>%
      ggplot() +
			ylim(0,NA) +
			xlab("Age (years)") +
			facet_wrap(~Sex,ncol=n_col)

      p = p + ylab("Length (cm)") +
              geom_hline(yintercept=0)
      if(input$growth_var)
      {
        p = p + geom_ribbon(aes(x=age,ymin=l95,ymax=u95,group=id,fill=id),alpha=0.25)
      }
	  p = p + geom_path(aes(x=age,y=pred_len,group=id,color=id),linewidth=1.5)
    
      p = p + viridis::scale_color_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			viridis::scale_fill_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
							panel.grid.major = element_line(color = 'gray70',linetype = "dotted"), 
							panel.grid.minor = element_line(color = 'gray70',linetype = "dotted"),
							strip.background =element_rect(fill="white"),
							legend.key = element_rect(fill = "white"))
			
    return(p)
  },
  height=function(){
        return(height_per_panel*1.5)
  })
# m_plot
  output$m_plot = renderPlot({
    input_models = unique(filtered_id())
    input_bio_sex = input$bio_sex
    if(length(input_models) < 1|length(input_bio_sex)<1){
      return()
    }
    selected_models = sapply(input_models,function(x)all_dirs[grep(x,all_dirs,fixed=TRUE)])
    
    plot_dt = rbindlist(lapply(selected_models,function(x)fread(file.path(model_stem,x,"bio.csv")))) %>%
             .[,.(id,Sex,M_age0,M_nages)] %>%
             melt(.,id.vars=c("id","Sex")) %>%
             setnames(.,c("variable","value"),c("age_cat","m"))

    max_age_dt = rbindlist(lapply(selected_models,function(x)fread(file.path(model_stem,x,"selex_a.csv"))[,.(M_age0=min(as.numeric(variable)),M_nages=max(as.numeric(variable))),by=id])) %>%
    melt(.,id.vars=c("id")) %>%
    setnames(.,c("variable","value"),c("age_cat","age"))


    plot_dt = merge(plot_dt,max_age_dt,by=c("id","age_cat"),allow.cartesian=TRUE) %>%
              .[,Sex:=factor(Sex,levels=as.character(2:1),labels=c("Male","Female"))] %>%
              .[Sex %in% input_bio_sex]

    if(nrow(plot_dt) == 0){
      return()
    }
    
    n_col = length(input_bio_sex)
    
    p = plot_dt %>%
      ggplot() +
			ylim(0,NA) +
			xlab("Age (years)") +
			facet_wrap(~Sex,ncol=n_col)

      p = p + ylab("Natural mortality (M)") +
              geom_hline(yintercept=0)
	  p = p + geom_path(aes(x=age,y=m,group=id,color=id),linewidth=1.5) +
            geom_point(aes(x=age,y=m,group=id,fill=id),shape=21,size=5)
      p = p + viridis::scale_color_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			viridis::scale_fill_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
							panel.grid.major = element_line(color = 'gray70',linetype = "dotted"), 
							panel.grid.minor = element_line(color = 'gray70',linetype = "dotted"),
							strip.background =element_rect(fill="white"),
							legend.key = element_rect(fill = "white"))
			
    return(p)
  },
  height=function(){
        return(height_per_panel*1.5)
  })
  # lw_plot
  output$lw_plot = renderPlot({
    input_models = unique(filtered_id())
    input_bio_sex = input$bio_sex
    if(length(input_models) < 1|length(input_bio_sex)<1){
      return()
    }
    selected_models = sapply(input_models,function(x)all_dirs[grep(x,all_dirs,fixed=TRUE)])
    
    plot_dt = rbindlist(lapply(selected_models,function(x)fread(file.path(model_stem,x,"bio.csv")))) %>%
             .[,.(id,Sex,WtLen1,WtLen2)] 

    max_len_dt = rbindlist(lapply(selected_models,function(x)fread(file.path(model_stem,x,"selex_l.csv"))[,.(max_len=max(as.numeric(variable))),by=id])) 
    len_dt.list = as.list(rep(NA,length(selected_models)))
    for(i in 1:length(selected_models)){
        len_dt.list[[i]] = as.data.table(expand.grid(id=input_models[i],len=0:ceiling(max_len_dt[id==input_models[i]]$max_len)))
    }
    len_dt = rbindlist(len_dt.list)

    plot_dt = merge(plot_dt,len_dt,by="id",allow.cartesian=TRUE) %>%
              .[,pred_weight := WtLen1*len^WtLen2] %>%
              .[,.(id,Sex,len,pred_weight)] %>%
              .[,Sex:=factor(Sex,levels=as.character(2:1),labels=c("Male","Female"))] %>%
              .[Sex %in% input_bio_sex]

    if(nrow(plot_dt) == 0){
      return()
    }
    
    n_col = length(input_bio_sex)
    
    p = plot_dt %>%
      ggplot() +
			ylim(0,NA) +
			xlab("Length (cm)") +
			facet_wrap(~Sex,ncol=n_col)

      p = p + ylab("Weight (kg)") +
              geom_hline(yintercept=0)
	  p = p + geom_path(aes(x=len,y=pred_weight,group=id,color=id),linewidth=1.5)
      p = p + viridis::scale_color_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			viridis::scale_fill_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
							panel.grid.major = element_line(color = 'gray70',linetype = "dotted"), 
							panel.grid.minor = element_line(color = 'gray70',linetype = "dotted"),
							strip.background =element_rect(fill="white"),
							legend.key = element_rect(fill = "white"))
			
    return(p)
  },
  height=function(){
        return(height_per_panel*1.5)
  })

  # catch_total_plot
  output$catch_total_plot = renderPlot({
    input_models = unique(filtered_id())
    input_units = input$catch_units
    input_aggregate = input$catch_by_fleet
    input_fill = input$catch_fill
    input_rescale = input$catch_rescale
    input_chart = input$catch_chart
    if(length(input_models) < 1){
      return()
    }
    selected_models = sapply(input_models,function(x)all_dirs[grep(x,all_dirs,fixed=TRUE)])
    # check to make sure catch.csv exists for all the selected models and update if needed
    keep_selected = which(sapply(file.path(model_stem,selected_models,"catch.csv"),file.exists)==TRUE)
    if(length(keep_selected)==0){
      return(warning("The file catch.csv does not exist for any of the files you selected."))
    } else {
      selected_models = selected_models[keep_selected]
    }
    plot_dt = rbindlist(lapply(selected_models,function(x)fread(file.path(model_stem,x,"catch.csv")))) %>%
             .[,.(id,Fleet,Time,Obs,Exp,sel_bio,sel_num)] %>%
             setnames(.,c("Fleet","Time","Obs","Exp"),c("fleet","time","obs","exp"))
    u_id = unique(plot_dt$id)
    # read in fleet_summary to get names and units
    fleet_dt.list = as.list(rep(NA,length(selected_models)))
    for(i in 1:length(selected_models)){
      fleet_dt.list[[i]] = fread(file.path(model_stem,selected_models[i],"fleet_summary.csv")) %>%
                           .[,id:=u_id[i]] %>%
                           .[,.(id,fleet,fleetname,type,units)]

    }
    fleet_dt = rbindlist(fleet_dt.list)
    plot_dt = merge(plot_dt,fleet_dt,by=c("id","fleet")) %>%
              .[,id2:=paste0(fleetname,"_",units)]

    if(nrow(plot_dt) == 0){
      return()
    }

    # plot if models == 1
    if(length(selected_models)==1){
      # define fill
      p = plot_dt[,.(id,fleet,fleetname,units,time,sel_bio,sel_num)] 
      if(input_fill == "None"){
        p = p %>% .[,fill:=id]
      } else if(input_fill == "Fleet"){
         p = p %>% .[,fill:=fleetname] %>% 
            .[,fill:=factor(fill,levels=unique(p$fleetname))]
      } else {
         p = p %>% .[,fill:=c("mt","n (1000s)")[units]]
      }
      # define catch
      if(input_aggregate == "TRUE"){
        if(input_units == "Numbers"){
          p = p[,.(catch=sum(sel_num)),by=.(time,fill)]
          ylab_txt = "Catch (1000s N)"
        } else {
          p = p[,.(catch=sum(sel_bio)),by=.(time,fill)]
          ylab_txt = "Catch (mt)"
        }
        if(input_rescale == "TRUE"){
          p = p %>% 
            .[,catch:=catch/mean(catch)]
          ylab_txt = paste0(ylab_txt,": re-scaled")
        }
      } else {
        if(input_units == "Numbers"){
          p = p %>% 
            .[,facet:=factor(fleetname,levels=unique(p$fleetname))] %>%
            .[,.(catch=sum(sel_num)),by=.(time,facet,fill)]
            ylab_txt = "Catch (1000s N)"
        } else {
          p = p %>% 
            .[,facet:=factor(fleetname,levels=unique(p$fleetname))] %>%
            .[,.(catch=sum(sel_bio)),by=.(time,facet,fill)]
          ylab_txt = "Catch (mt)"
        }
        if(input_rescale == "TRUE"){
          p = p %>% 
            .[,catch:=catch/mean(catch),by=facet]
          ylab_txt = paste0(ylab_txt,": re-scaled")
        }
      }
      p = p %>%
        ggplot() +
			  ylim(0,NA) +
			  xlab("Time")
        if(input_aggregate=="FALSE"){
          n_row = 5
          p = p + facet_wrap(~facet,nrow=n_row)
        }
      p = p +
          ylab(ylab_txt)
        if(input_chart == "Bar"){
          p = p + geom_bar(aes(x=time,y=catch,fill=fill),position="stack",stat = "identity",color="white")
        } else {
          p = p + geom_area(aes(x=time,y=catch,fill=fill))
        }
      p = p + geom_hline(yintercept=0) + 
          viridis::scale_color_viridis("Category",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			    viridis::scale_fill_viridis("Category",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			    theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
							panel.grid.major = element_line(color = 'gray70',linetype = "dotted"), 
							panel.grid.minor = element_line(color = 'gray70',linetype = "dotted"),
							strip.background =element_rect(fill="white"),
							legend.key = element_rect(fill = "white"))
      return(p)
    } else {
        if(length(selected_models)>1&input_aggregate=="FALSE"&uniqueN(plot_dt$id2) != max(plot_dt$fleet)){
            return()
        }
        p = plot_dt[,.(id,fleet,fleetname,units,time,sel_bio,sel_num)] 
        p = p %>%
            .[,id:=factor(id,levels=unique(p$id))]
        if(input_units == "Numbers"){
          p = p[,.(catch=sum(sel_num)),by=.(id,time)]
          ylab_txt = "Catch (1000s N)"
        } else {
          p = p[,.(catch=sum(sel_bio)),by=.(id,time)]
          ylab_txt = "Catch (mt)"
        }
        if(input_rescale == "TRUE"){
          p = p %>% 
            .[,catch:=catch/mean(catch),by=.(id)]
          ylab_txt = paste0(ylab_txt,": re-scaled")
        }

         p = p %>% .[order(id,time)] %>%
        ggplot() +
			  ylim(0,NA) +
			  xlab("Time")
      p = p +
          ylab(ylab_txt)
          p = p + geom_path(aes(x=time,y=catch,color=id),linewidth=2)
        
      p = p + geom_hline(yintercept=0) + 
          viridis::scale_color_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			    viridis::scale_fill_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
			    theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
							panel.grid.major = element_line(color = 'gray70',linetype = "dotted"), 
							panel.grid.minor = element_line(color = 'gray70',linetype = "dotted"),
							strip.background =element_rect(fill="white"),
							legend.key = element_rect(fill = "white"))
      return(p)
    }
  },
  height=function(){
    if(input$catch_by_fleet)
    {
      return((height_per_panel*1.5))
    } else {
      return((height_per_panel*1.5))
    }
  })


  # catch_total_plot
  output$cpue_obs = renderPlot({
    input_models = unique(filtered_id())
    input_se = input$cpue_se
    input_multi = input$cpue_multi
    input_fit = input$cpue_fit
    input_log = input$cpue_log
    input_varadj = input$cpue_varadj
    input_lambda = input$cpue_lambda

    if(length(input_models) < 1){
      return()
    }
    selected_models = sapply(input_models,function(x)all_dirs[grep(x,all_dirs,fixed=TRUE)])
    # check to make sure cpue.csv exists for all the selected models and update if needed
    keep_selected = which(sapply(file.path(model_stem,selected_models,"cpue.csv"),file.exists)==TRUE)
    if(length(keep_selected)==0){
      return(warning("The file cpue.csv does not exist for any of the files you selected."))
    } else {
      selected_models = selected_models[keep_selected]
    }
    plot_dt = rbindlist(lapply(selected_models,function(x)fread(file.path(model_stem,x,"cpue.csv")))) %>%
             .[,.(id,Fleet,Fleet_name,Time,Obs,Exp,SE)] %>%
             setnames(.,c("Fleet","Fleet_name","Time","Obs","Exp","SE"),c("fleet","fleetname","time","obs","exp","se")) %>%
             .[,obs:=round(obs,digits=3)] %>%
             .[,id2:=as.numeric(as.factor(paste0(fleetname,"_",obs)))] %>%
             .[,id3:=as.numeric(as.factor(paste0(id,"_",fleetname)))]
    
    u_id = unique(plot_dt$id)
    # read in varadj and lambdas
    varlambda_dt.list = as.list(rep(NA,length(selected_models)))
    for(i in 1:length(selected_models)){
      varlambda_dt.list[[i]] = 

      tmp_var = fread(file.path(model_stem,selected_models[i],"varadj.csv"))
      if(nrow(tmp_var)==0){
        tmp_var = data.table(fleet=unique(plot_dt[id==u_id[i]]$fleet)) %>%
                  .[,id:=u_id[i]] %>%
                  .[,factor:=1] %>%
                  .[,value:=0] %>%
                  .[,.(id,factor,fleet,value)]
      } else {
        tmp_var = tmp_var[factor==1]
        tmp_var2 = data.table(fleet=unique(plot_dt[id==u_id[i]]$fleet)) %>%
                  .[,id:=u_id[i]] %>%
                  .[,factor:=1] %>%
                  .[,value:=0] %>%
                  .[,.(id,factor,fleet,value)] %>%
                  .[!(fleet%in%tmp_var$fleet)]
        tmp_var = rbind(tmp_var,tmp_var2) %>% .[order(fleet)]
      }
      tmp_var = tmp_var %>%
                .[,.(id,fleet,value)] %>%
                setnames(.,c("fleet","value"),c("fleet","varadj"))

      tmp_lambda = fread(file.path(model_stem,selected_models[i],"lambdas.csv")) %>%
                  .[like_comp==1] %>%
                  .[,.(id,fleet,value)]
      if(nrow(tmp_var)==0){
        tmp_lambda = data.table(fleet=unique(plot_dt[id==u_id[i]]$fleet)) %>%
                  .[,id:=u_id[i]] %>%
                  .[,value:=0] %>%
                  .[,.(id,fleet,value)]
      }  else {
        tmp_lambda2 = data.table(fleet=unique(plot_dt[id==u_id[i]]$fleet)) %>%
                  .[,id:=u_id[i]] %>%
                  .[,value:=0] %>%
                  .[,.(id,fleet,value)] %>%
                  .[!(fleet%in%tmp_lambda$fleet)]
        tmp_lambda = rbind(tmp_lambda,tmp_lambda2) %>% .[order(fleet)]
      }
      tmp_lambda = tmp_lambda %>%
                setnames(.,c("value"),c("lambda"))        

      varlambda_dt.list[[i]] = merge(tmp_var,tmp_lambda,by=c("id","fleet"))
    }
    varlambda_dt = rbindlist(varlambda_dt.list)
    plot_dt = merge(plot_dt,varlambda_dt,by=c("id","fleet")) 

    if(input_varadj == "FALSE"){
      plot_dt = plot_dt %>%
                .[,se:=se-varadj]
    }

    if(input_log == "TRUE"){
      plot_dt = plot_dt %>%
                .[,obs:=log(obs)] %>%
                .[,exp:=log(exp)] %>%
                .[,lse:=obs-se] %>%
                .[,use:=obs+se]
                ylab_txt = "Index (log-scale)"
                yint = 0
    } else {
      plot_dt = plot_dt %>%
                .[,lse:=exp(log(obs)-se)] %>%
                .[,use:=exp(log(obs)+se)]
                ylab_txt = "Index"
                yint = 1
    }

    if(input_lambda == "TRUE"){
      plot_dt = plot_dt[lambda==1]
    }

    if(nrow(plot_dt) == 0){
      return()
    }

    # plot if models == 1
    if(length(selected_models)==1){
      p = plot_dt %>% .[order(id3,time)] %>%
        ggplot()
        if(input_log != "TRUE"){
          p = p + ylim(0,NA)
        }
			  p = p + xlab("Time") +
        ylab(ylab_txt) +
        geom_hline(yintercept=yint,linetype="dashed",linewidth=1) 
      if(input_multi == "TRUE"){
        p = p + facet_wrap(~fleetname)
        if(input_se == "TRUE"){
        p = p + geom_segment(aes(x=time,xend=time,y=lse,yend=use),linewidth=1)
        }
        p = p + geom_point(aes(x=time,y=obs),fill="gray50",shape=21,size=3)
        if(input_fit == "TRUE"){
          p = p + geom_path(aes(x=time,y=exp,color=id),linewidth=1)
        }  
        p = p +  
            viridis::scale_color_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
            viridis::scale_fill_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
            theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
                panel.grid.major = element_line(color = 'gray70',linetype = "dotted"), 
                panel.grid.minor = element_line(color = 'gray70',linetype = "dotted"),
                strip.background =element_rect(fill="white"),
                legend.key = element_rect(fill = "white"))
      } else {
        p = p + geom_path(aes(x=time,y=obs,group=id3,color=fleetname),linewidth=1)
        p = p +  
            viridis::scale_color_viridis("Index",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
            viridis::scale_fill_viridis("Index",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
            theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
                panel.grid.major = element_line(color = 'gray70',linetype = "dotted"), 
                panel.grid.minor = element_line(color = 'gray70',linetype = "dotted"),
                strip.background =element_rect(fill="white"),
                legend.key = element_rect(fill = "white"))
      }
      return(p)
    } else {
        if(length(selected_models)>1&input_multi=="TRUE"& mean(table(plot_dt$id2) %% uniqueN(plot_dt$id) == 0)!= 1){
            return()
        }
        p = plot_dt %>% .[order(id3,time)] %>%
        ggplot()
        if(input_log != "TRUE"){
          p = p + ylim(0,NA)
        }
			  p = p + xlab("Time") +
        ylab(ylab_txt) +
        geom_hline(yintercept=yint,linetype="dashed",linewidth=1) 
      if(input_multi == "TRUE"){
        p = p + facet_wrap(~fleetname)
        if(input_se == "TRUE"){
        p = p + geom_segment(aes(x=time,xend=time,y=lse,yend=use,group=id3,color=id),linewidth=1)
        }
        p = p + geom_point(aes(x=time,y=obs),fill="gray50",shape=21,size=3)
        if(input_fit == "TRUE"){
          p = p + geom_path(aes(x=time,y=exp,color=id),linewidth=1)
        }  
        p = p +  
            viridis::scale_color_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
            viridis::scale_fill_viridis("Model",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
            theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
                panel.grid.major = element_line(color = 'gray70',linetype = "dotted"), 
                panel.grid.minor = element_line(color = 'gray70',linetype = "dotted"),
                strip.background =element_rect(fill="white"),
                legend.key = element_rect(fill = "white"))
      } else {
        p = p + geom_path(aes(x=time,y=obs,group=id3,color=fleetname),linewidth=1)
        p = p +  
            viridis::scale_color_viridis("Index",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
            viridis::scale_fill_viridis("Index",begin = 0.1,end = 0.8,direction = 1,option = "H",discrete=TRUE) +
            theme(panel.background = element_rect(fill = "white", color = "black", linetype = "solid"),
                panel.grid.major = element_line(color = 'gray70',linetype = "dotted"), 
                panel.grid.minor = element_line(color = 'gray70',linetype = "dotted"),
                strip.background =element_rect(fill="white"),
                legend.key = element_rect(fill = "white"))
      }
      return(p)
    }
  },
  height=function(){

      return((height_per_panel*1.5))

  })
} # End of server
