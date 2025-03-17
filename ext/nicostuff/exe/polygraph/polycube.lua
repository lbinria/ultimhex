-- Lua (keep this comment, it is an indication for editor's 'run' command)

-- TODO
-- periodicity: compute isometry, match vertices and use it in smoother
-- find paradigm for data management i.e. coherent data state between saved/on_disk/in_memory
-- quantization decaler tout les iso suivants en cas d'egalité
-- remove hex selection
-- interactive tri segmentation into charts
-- check local configuration of flagging
-- check validity of faces to pad



function set_gui_visibility(b)
  gom.set_environment_value('gui:presentation_mode','true')
  graphite_gui.visible = b
  Stats.visible= b
  scene_graph_gui.visible = b
  toolbox_gui.visible = b
  console_gui.visible = b
  text_editor_gui.visible =  b
  object_properties_gui.visible = b
  camera_gui.visible = b
  console_gui.visible = b
end


cull_front = false
prefered_paint_tool = "OGF::MeshGrobPaint" 
autoconfirm = false
use_gmsh_regions = false

auto_smooth = true
safe_mode = false

pickflag = false

function switch_interaction_mode()
  if state=="embedit" then
    if scene_graph.current_object == "quadchart" then 
      set_embedit_mode_picking()
    else 
      set_embedit_mode_painting()
    end
  end
if state =="paint" then
  if pickflag then 
    tools.tool(prefered_paint_tool)  
    scene_graph.current().shader.mesh_style = 'true; 1 1 1  1; 1'
  else 
    tools.tool('OGF::MeshGrobProbe')
    scene_graph.current().shader.mesh_style = 'false; 1 1 1  1; 1'
  end
  pickflag = not pickflag
end 
end


key_pressed = {}
function key_down(k)
  if k=="a" then switch_interaction_mode() end
    key_pressed[k] = true
end

function key_up(k)
    key_pressed[k] = false
end
gom.connect(main.render_area.key_down, key_down)
gom.connect(main.render_area.key_up, key_up)






bin_path = "C:/NICO/prog/nicostuff/build/RelWithDebInfo/"
path,project_name  = "",""
state = "idle"

function load_project()
  dos_path = path .. project_name
  dos_path = dos_path:gsub("/", "\\")
  files = {"tri","flag","hex"}
  for trash,name in pairs(files)  do
   os.execute("copy \"" .. dos_path .. "\\" .. name .. "_tmp.geogram\" \"" .. path .. project_name .. "/" .. name .. ".geogram\"" )
  end
  set_state("viewall")
end

function save_project()
  to_disk()
  dos_path = path .. project_name
  dos_path = dos_path:gsub("/", "\\")
  files = {"tri","flag","hex"}
  for trash,name in pairs(files)  do
   os.execute("copy \"" .. dos_path .. "\\" .. name .. ".geogram\" \"" .. path .. project_name .. "/" .. name .. "_tmp.geogram\"" )
  end
end


function filename(name)  return path .. project_name .. "/" .. name .. ".geogram" end

function file_exists(name)
   local f=io.open(path .. project_name .. "/" .. name .. ".geogram","r")
   if f~=nil then io.close(f) return true else return false end
end

function file_remove(name)
  cmd = "del  \"" .. path .. project_name .. "/" .. name .. ".geogram\""
  os.execute(cmd:gsub("/", "\\"))
end

function refresh_current(name)
  scene_graph.clear()
  scene_graph.load_object({value=filename(name),type='default',invoked_from_gui=false})
end

function run(binary_name)
  file_remove("feedback")
  os.execute(bin_path .. binary_name .. " projectpath=\"" .. path .. project_name .."\" run_from=graphite")
end

function set_current(name)
  if scene_graph.current_object ~= name then
    refresh_current(name)
    end
end
function to_disk()
  for i=0,scene_graph.nb_children-1 do
    name = scene_graph.ith_child(i).name
    print(name)
    scene_graph.current_object = name
    scene_graph.save_current_object(filename(name))
  end
end


function set_embedit_mode_picking()
  scene_graph.current_object = "quadchart"
  scene_graph.current().visible = false
  scene_graph.current_object ="trichart"
  if cull_front then 
    scene_graph.current().shader.culling_mode = "CULL_FRONT" 
  else 
    scene_graph.current().shader.culling_mode = "NO_CULL" 
  end
  scene_graph.current().visible = true
  tools.tool('OGF::MeshGrobProbe')
  scene_graph.current().shader.mesh_style = 'false; 1 1 1  1; 1'

end

function set_embedit_mode_painting()
  scene_graph.current_object = "trichart"
  scene_graph.current().visible = false
  scene_graph.current_object ="quadchart"
  scene_graph.current().visible = true
  if cull_front then 
    scene_graph.current().shader.culling_mode = "CULL_FRONT" 
  else 
    scene_graph.current().shader.culling_mode = "NO_CULL" 
  end
  tools.tool(prefered_paint_tool)
  scene_graph.current().shader.mesh_style = 'true; 1 1 1  1; 1'
end


function set_embedit_mode()
  scene_graph.load_object({value=filename("trichart"),type='default',invoked_from_gui=false})
  scene_graph.current().shader.painting= 'ATTRIBUTE'
  scene_graph.current().shader.attribute = 'facets.chart'
  scene_graph.current().shader.autorange()
  scene_graph.current().shader.mesh_style = 'true; 1 1 1  1; 1'
  scene_graph.current().shader.vertices_style = 'false; 0 1 0 1; 5'

  scene_graph.load_object({value=filename("quadchart"),type='default',invoked_from_gui=false})
  scene_graph.current().shader.painting= 'ATTRIBUTE'
  scene_graph.current().shader.attribute = 'facets.chart'
  scene_graph.current().shader.autorange()
  scene_graph.current().shader.mesh_style = 'true; 1 1 1  1; 1'
  scene_graph.current().shader.vertices_style = 'false; 0 1 0 1; 5'
  set_embedit_mode_picking()
end

function set_viewall_mode()
  files = {"hex","flag","tri"}
  for trash,name in pairs(files)  do
      if (file_exists(name)) then
         if name ~="hex" then
          scene_graph.load_object({value=filename(name),type='default',invoked_from_gui=false})
         else
          run("poly_hexedit.exe algo=hexview ")
          scene_graph.load_object({value=filename("singuview"),type='default',invoked_from_gui=false})
          scene_graph.current().shader.painting= 'ATTRIBUTE'
          scene_graph.current().shader.edges_style = 'true; 0 0 0  1; 10'
          scene_graph.current().shader.vertices_style = 'false; 0 1 0 1; 5'
          scene_graph.current().shader.attribute = 'edges.singu'
          scene_graph.current().shader.attribute_min = 2
          scene_graph.current().shader.attribute_max = 4
          scene_graph.current().visible = false
          scene_graph.load_object({value=filename("hexview"),type='default',invoked_from_gui=false})
          scene_graph.current().shader.painting= 'ATTRIBUTE'
          scene_graph.current().shader.attribute = 'cells.block'
          scene_graph.current().shader.autorange()
          scene_graph.current().shader.edges_style = 'true; 0.5 0.5 1  1; 1'
          scene_graph.current().shader.shrink = 1
          scene_graph.current().shader.vertices_style = 'false; 0 1 0 1; 5'
         end
         return
        end
      end
end


function set_state(new_state)
  if state == new_state then return end
  tools.tool('OGF::GrobPan')
  state = new_state
  scene_graph.clear()
  if state== "viewall" then
    set_viewall_mode()
  end

  print(state)
  if state == "paint" then
     refresh_current("flag")
     set_paint_axis(current_paint_value)
  end

  if state == "hexcollapse" then
     refresh_current("hex")
     set_layer_collapse_mode()
  end

  if state == "hexpad" then
    run("poly_hexedit.exe algo=padInterface ")
    set_pad_mode()
  end

  if state == "embedit" then
        if use_gmsh_regions then run("poly_hexedit.exe algo=embeditinit gmshchart=true") 
        else run("poly_hexedit.exe algo=embeditinit  gmshchart=false") end
    set_embedit_mode()
  end

  if state == "killhex" then
    run("poly_hexedit.exe algo=killhexinit ")
     set_current("hex")
    scene_graph.current().shader.mesh_style = 'true; 0 0 0 1; 1'
  scene_graph.current().shader.painting = 'ATTRIBUTE'
  scene_graph.current().shader.attribute = 'cells.tokill'
  scene_graph.current().shader.attribute_min = 0
  scene_graph.current().shader.attribute_max = 1
 tools.tool(prefered_paint_tool)
 tools.current().autorange = false
 tools.current().accumulate = false
  tools.current().value = 1

  end
  

  if state =="paintfeedback" then
    set_current("flag")
    scene_graph.load_object({value=filename("feedback"),type='default',invoked_from_gui=false})
    scene_graph.current().shader.painting= 'ATTRIBUTE'
    scene_graph.current().shader.attribute = 'edges.fail_id'
    scene_graph.current().shader.attribute_min = -1
    scene_graph.current().shader.attribute_max = 0
  
    scene_graph.current().shader.edges_style = 'true; 1 0 .5  1; 20'
    scene_graph.current().shader.vertices_style = 'false; 1 0 0 1; 0'
    file_remove("feedback")
  end



  -- if state =="hexcollapsefeedback" then
  --   scene_graph.set_current("hex")
  --   scene_graph.current().shader.shrink = 5
  --   scene_graph.current().shader.vertices_style = 'false; 0 1 0 1; 5'
  --   print("collapse failed, see problematic facets in feedback")
  --   scene_graph.load_object({value=filename("feedback"),type='default',invoked_from_gui=false})
  --   file_remove("feedback")
  -- end

end
function forced_set_state(new_state)
  state = "idle"
  set_state(new_state)
end



--                      _
--                     | |
--   ___ _ __ ___  __ _| |_ ___
--  / __| '__/ _ \x/ _` | __/ _ \
-- | (__| | |  __/ (_| | ||  __/
--  \___|_|  \___|\__,_|\__\___|

function header_gui()
  if imgui.Button('create/load project from file') then
   tools.tool('OGF::GrobPan')
   local file = scene_graph.current().filename
   path = file:match("(.*/)")
   project_name = string.match(file,"/([^/]*)[.]")
   if project_name ~= "tri"  and  project_name ~= "flag" and project_name ~= "hex" and project_name ~= "quad" then
    os.execute("mkdir \"" .. path .. project_name .. "\"")
    scene_graph.save_current_object(filename("tri"))
    load_project()
   end
   set_state("viewall")
 end
end



function bottom_gui()
imgui.PushStyleColor(ImGuiCol_Button,0xFFFFFFFF)

  if imgui.Button('Save') then
    save_project()
  end
  imgui.SameLine()
  if imgui.Button('(Re)Load') then
    load_project()
  end
  if imgui.Button("scenegraph2disk") then
    to_disk()
  end
 imgui.PopStyleColor(1)

end


--              _       _
--             (_)     | |
-- _ __    __ _  _ _ __| |_
-- | '_ \ / _` | | '_ \| __|
-- | |_) | (_| | | | | | |_
-- | .__/ \__,_|_|_| |_|\__|
-- | |
-- |_|

flag_step = 5
const_decal=2
current_paint_value = -1
function set_paint_axis(axis)
  scene_graph.current().shader.mesh_style = 'true; 0 0 0 1; 1'
  scene_graph.current().shader.lighting=false
  scene_graph.current().shader.painting = 'ATTRIBUTE'
  scene_graph.current().shader.attribute = 'facets.flag'
  scene_graph.current().shader.attribute_min = -1
  scene_graph.current().shader.attribute_max = flag_step*5+const_decal
 current_paint_value = axis
 tools.tool(prefered_paint_tool)
 tools.current().autorange = false
 tools.current().accumulate = false
  tools.current().value = -1
 if axis~=-1 then
  tools.current().value = axis*flag_step + const_decal
 end
end




function auto_fill_flag()
  tools.tool('OGF::GrobPan')
  to_disk()
  run("poly_flag.exe ALGO=fill")
  state = "idle"
  set_state("paint")
end

nhexwanted = 1000



function paintfeedback_gui()
    imgui.Text("Degenerated flag configurations\nare shown in 3D windows.\nIs it fine for you ?")
    if file_exists("hex") then
      imgui.Text("Hex mesh already exists.\nThis opration will erase it.")
    end
    if imgui.Button('Cancel##hexgen') then set_state('paint') end
      imgui.SameLine()
      if imgui.Button('Confirm##hexgen') or autoconfirm then
        tools.tool('OGF::GrobPan')
        run("poly_cube.exe nbhex=" .. nhexwanted)
        set_state( "hexcollapse")
        if use_gmsh_regions then run("poly_hexedit.exe algo=embeditinit gmshchart=true") 
        else run("poly_hexedit.exe algo=embeditinit  gmshchart=false") end
  end
end


function paint_gui()


  imgui.Text("Tools")
  imgui.PushStyleColor(ImGuiCol_Button,0xFFC55e17)
  if imgui.Button('paint-X') then set_paint_axis(0) end imgui.SameLine()
  imgui.PushStyleColor(ImGuiCol_Button,0xFF789608)
  if imgui.Button('paint+X') then set_paint_axis(1) end
  imgui.PushStyleColor(ImGuiCol_Button,0xFF2ef877)
  if imgui.Button('paint-Y') then set_paint_axis(2) end imgui.SameLine()
  imgui.PushStyleColor(ImGuiCol_Button,0xFF2ef8bb)
  if imgui.Button('paint+Y') then set_paint_axis(3) end
  imgui.PushStyleColor(ImGuiCol_Button,0xFF1ee1f9)
  if imgui.Button('paint-Z') then set_paint_axis(4) end imgui.SameLine()
  imgui.PushStyleColor(ImGuiCol_Button,0xFF1e56f9)
  if imgui.Button('paint+Z') then set_paint_axis(5) end
  imgui.PushStyleColor(ImGuiCol_Button,0xFFfa1d03)
  if imgui.Button('ERASE CONSTRAINTS') then set_paint_axis(-1) end
  imgui.PopStyleColor(7)

 imgui.Separator()
 imgui.Text("Commands")


 if imgui.Button('BEND FROM CONSTRAINT') then  auto_fill_flag() end
 if imgui.Button('LOCK BEND') then  
  to_disk()
  run("poly_flag.exe ALGO=lockbend")
  state = "idle"
  set_state("paint")
 end
 if imgui.Button('UNLOCK BEND') then  
  to_disk()
  run("poly_flag.exe ALGO=unlockbend")
  state = "idle"
  set_state("paint")
 end



 imgui.Separator()
 imgui.Text("produce hex mesh")
  local sel,new_value = imgui.InputInt("nhex",nhexwanted)
  nhexwanted=new_value
 if imgui.Button('GENERATE HEX') then
    to_disk()
    run("poly_flag.exe ALGO=diagnostic")
    set_state("paintfeedback")
  end
end


--  _                          _ _ _
-- | |                        | (_) |
-- | |__   _____  __   ___  __| |_| |_
-- | '_ \ / _ \ \/ /  / _ \/ _` | | __|
-- | | | |  __/>  <  |  __/ (_| | | |_
-- |_| |_|\___/_/\_\  \___|\__,_|_|\__|

hex_tool = ""

function set_default_hex_shader()
  tools.tool('OGF::GrobPan')
 scene_graph.current().shader.lighting=true
 scene_graph.current().shader.painting = 'SOLID_COLOR'
  scene_graph.current().shader.mesh_style = 'true; 0 0 0 1; 1'
  scene_graph.current().shader.vertices_style = 'false; 0 1 0 1; 5'
end


function set_layer_collapse_mode()
  scene_graph.current().shader.shrink = 0
 scene_graph.current().shader.lighting=false
 scene_graph.current().shader.painting = 'SOLID_COLOR'
  scene_graph.current().shader.lighting = false
  scene_graph.current().shader.mesh_style = 'true; 0 0 0 1; 1'
  scene_graph.current().shader.edges_style = 'true; 0 0 1 1; 5'
  scene_graph.current().shader.vertices_style = 'true; 0 1 0 1; 5'
  tools.tool('OGF::MeshGrobCreateEdge')
end


pad_mode_inside = false
function set_pad_mode()
  to_disk()
  if pad_mode_inside then
   set_current("padselectin")
   scene_graph.current().shader.culling_mode= 'CULL_FRONT'
  else
   set_current("padselecton")
  end
  scene_graph.current().shader.vertices_style = 'false; 0 1 0 1; 5'
  scene_graph.current().shader.painting = 'ATTRIBUTE'
  scene_graph.current().shader.attribute = 'facets.fpaint'
  scene_graph.current().shader.mesh_style = 'true; 0 0 0 1; 1'
  scene_graph.current().shader.attribute_min = 0
  scene_graph.current().shader.attribute_max = 1

  tools.tool(prefered_paint_tool)
  tools.current().autorange = false
  tools.current().accumulate = false
  tools.current().value = 1
end


 npadlayers =1




  ------------------------------------------------------------------------------------------------
  function hexcollapse_gui()

   if imgui.Button('Clear edges to collapse') then
     run("poly_hexedit.exe algo=clearlayerkill ")
     forced_set_state("hexcollapse")
    end

   if imgui.Button('Collapse') then
     to_disk()
     run("poly_hexedit.exe algo=layerkill ")
     forced_set_state("hexcollapse")

     if auto_smooth then smooth() end 

    end
--    imgui.Separator()
 --   if imgui.Button('Smooth') then
 --  to_disk()
--   os.execute(bin_path .. "poly_hexedit.exe algo=smooth projectpath=\"" .. path .. project_name .."\" run_from=graphite")
--   forced_set_state("hexcollapse")
--  end

   end



   pad_invert = false

   function hexpad_gui()
   if imgui.Button('edit on') then
     pad_mode_inside =false
     set_pad_mode()
   end
   imgui.SameLine()
   if imgui.Button('edit inside') then
     pad_mode_inside =true
     set_pad_mode()
   end

  
  local sel,new_val =imgui.Checkbox("InvertNormal##embedit",pad_invert)
  pad_invert = new_val
  local sel,new_value = imgui.InputInt("nPadLayers",npadlayers)
  npadlayers=new_value

   if imgui.Button('view pad faces') then
     to_disk()
     if pad_invert then  run("poly_hexedit.exe algo=padfaces invertNormal=true")
     else run("poly_hexedit.exe algo=padfaces invertNormal=false") end
     tools.tool('OGF::GrobPan')
     set_current("hex")
     scene_graph.current().shader.painting= 'SOLID_COLOR'
     scene_graph.current().shader.shrink = 8
     scene_graph.current().shader.vertices_style = 'false; 0 1 0 1; 5'
     scene_graph.current().shader.volume_style = 'true; 1 1 1 1'
     scene_graph.load_object({value=filename("feedback"),type='default',invoked_from_gui=false})
     scene_graph.current().shader.painting= 'SOLID_COLOR'
     scene_graph.current().shader.surface_style = " true ; 0 0 1 1"
     scene_graph.current().shader.two_sided= 'true'

     scene_graph.load_object({value=filename("feedbackerror"),type='default',invoked_from_gui=false})
     scene_graph.current().shader.painting= 'SOLID_COLOR'
     scene_graph.current().shader.vertices_style = 'false; 0 0 1 1; 5'
     scene_graph.current().shader.edges_style = 'true; 1 0 0 1; 30'
     
     scene_graph.load_object({value=filename("cuts"),type='default',invoked_from_gui=false})
     scene_graph.current().shader.painting= 'ATTRIBUTE'
     scene_graph.current().shader.attribute = "edges.out_attr"
      scene_graph.current().shader.attribute_min = 1
    scene_graph.current().shader.attribute_max = 4
     scene_graph.current().shader.vertices_style = 'false; 0 0 1 1; 5'
     scene_graph.current().shader.edges_style = 'true; 0 0 1 1; 2'
     
   end



   if imgui.Button('apply') then
     to_disk()
     if pad_invert then  run("poly_hexedit.exe algo=pad nbpadlayers=" .. npadlayers .. " invertNormal=true")
     else run("poly_hexedit.exe algo=pad nbpadlayers=" .. npadlayers .. " invertNormal=false") end
     forced_set_state("hexpad")
     if auto_smooth then smooth() end
        
     end

   imgui.SameLine()
   if imgui.Button('Clear') then
    forced_set_state("hexpad")
  end
    -- if imgui.Button("do something") then
    --   warningmsg = "do you really wanna quit ?!?"
    --   print(cmd_confirmed)
    -- end
  imgui.Separator()
    if imgui.Button('Smooth') then
   to_disk()
   os.execute(bin_path .. "poly_hexedit.exe algo=smooth projectpath=\"" .. path .. project_name .."\" run_from=graphite")
   forced_set_state("hexpad")
  end

  end





--                _              _ _ _
--               | |            | (_) |
--  ___ _ __ ___ | |__   ___  __| |_| |_
-- / _ \ '_ ` _ \| '_ \ / _ \/ _` | | __|
--|  __/ | | | | | |_) |  __/ (_| | | |_
-- \___|_| |_| |_|_.__/ \___|\__,_|_|\__|


function embedit_gui()

imgui.Text("Use F1 to switch between\n Pick/Paint tools...\n if properly configured :)")
  if imgui.Button("Pick##embedit") then set_embedit_mode_picking() end
  imgui.SameLine()
  if imgui.Button("Paint##embedit") then set_embedit_mode_painting() end
  if imgui.Button("Apply##embedit") then
    tools.tool('OGF::GrobPan')
    to_disk()
    run("poly_hexedit.exe algo=embeditapply ")
  if auto_smooth then smooth() end
  set_embedit_mode_painting()
end
  
---  if imgui.Button("smooth##embedit") then
--    tools.tool('OGF::GrobPan')
 --   to_disk()
--    os.execute(bin_path .. "poly_hexedit.exe algo=smooth projectpath=\"" .. path .. project_name .."\" run_from=graphite")
--    forced_set_state("embedit")
--   end
end


function killhex_gui()
  if imgui.Button("Apply##killhex") then
    to_disk()
    run("poly_hexedit.exe algo=killhexapply ")
     forced_set_state("killhex")

     if auto_smooth then smooth() end
  end
end

--        _                    _ _
--       (_)                  | | |
-- __   ___  _____      ____ _| | |
-- \ \ / / |/ _ \ \ /\ / / _` | | |
--  \ V /| |  __/\ V  V / (_| | | |
--   \_/ |_|\___| \_/\_/ \__,_|_|_|


viewall_gui_wait_confirmation = false

function viewall_gui()
  if viewall_gui_wait_confirmation then
    imgui.Text("WARNING :\nThis operation will reset\nthe flaggin.\nPrevious flagging will be lost.")
    if imgui.Button('cancel') then viewall_gui_wait_confirmation = false end
    imgui.SameLine()
    if imgui.Button('confirm') or autoconfirm then
      run("poly_flag.exe ALGO=unlockbend")
      auto_fill_flag()
      viewall_gui_wait_confirmation = false
     end
    return
  end
  if imgui.Button('init flag') then
    if file_exists("flag") then
        viewall_gui_wait_confirmation = true
    else
      auto_fill_flag()

  end
  end
--  if imgui.Button("smooth") then
--   to_disk()
--   os.execute(bin_path .. "poly_hexedit.exe algo=smooth projectpath=\"" .. path .. project_name .."\" run_from=graphite")
--   forced_set_state("viewall")
--  end
end



--                  _
--                 (_)
--  _ __ ___   __ _ _ _ __
-- | '_ ` _ \ / _` | | '_ \
-- | | | | | | (_| | | | | |
-- |_| |_| |_|\__,_|_|_| |_|

my_dialog = {}
my_dialog.visible = true
my_dialog.name = 'Polycube'
my_dialog.x ,my_dialog.y = 100,400

polycube_option_dialog = {}
polycube_option_dialog.visible = true
polycube_option_dialog.name = 'PolycubeOption'
polycube_option_dialog.x ,polycube_option_dialog.y = 200,600




function my_dialog.draw_window()
  
  --print(key_pressed['shift'])
  

  imgui.Text('POLYCUBE')
  imgui.Spacing()
  imgui.Separator()

 if project_name=="" then  header_gui()  return end



 if state == "paintfeedback" then  paintfeedback_gui() return end


 accessible_states = "viewall;paint;hexcollapse;hexpad;embedit;killhex"
 if not file_exists("hex") then accessible_states = "viewall;paint" end
 if not file_exists("flag") then accessible_states = "viewall" end
 local sel,new_state = imgui.Combo("STEP",state,accessible_states)

 imgui.Separator()
 set_state(new_state)
 if state == "viewall" then viewall_gui() end
 if state == "paint" then paint_gui() end
 if state == "hexcollapse" then  hexcollapse_gui() end
 if state == "hexpad" then  hexpad_gui() end
 if state == "embedit" then  embedit_gui() end
 if state == "killhex" then  killhex_gui() end
 
end

function smooth()
   to_disk()
   os.execute(bin_path .. "poly_hexedit.exe algo=smooth projectpath=\"" .. path .. project_name .."\" run_from=graphite")
   forced_set_state(state)
end

function polycube_option_dialog.draw_window()
  imgui.Text("Smoothing")

  local sel,new_val =imgui.Checkbox("auto_smooth",auto_smooth)
  auto_smooth = new_val

  if imgui.Button("smooth") then
  smooth()
  end


  imgui.Separator() 
  imgui.Text("OPTIONS")
  if scene_graph.current() ~= nil then

  local sel,new_val =imgui.Checkbox("CULL##embedit",cull_front)
  cull_front = new_val
  
  if cull_front then scene_graph.current().shader.culling_mode = "CULL_FRONT" 
  else scene_graph.current().shader.culling_mode = "NO_CULL" end

  local sel,new_val =imgui.Checkbox("use_gmsh_regions",use_gmsh_regions)
  use_gmsh_regions = new_val


  local sel,new_paint_tool = imgui.Combo("prefered paint tool",prefered_paint_tool,"OGF::MeshGrobPaint;OGF::MeshGrobPaintRect;OGF::MeshGrobPaintFreeform")
  prefered_paint_tool = new_paint_tool

  local sel,new_val =imgui.Checkbox("autoconfirm",autoconfirm)
  autoconfirm = new_val

  local sel,new_val =imgui.Checkbox("SAFE MODE",safe_mode)
  if safe_mode ~=new_val then set_gui_visibility(safe_mode) end
  safe_mode = new_val


  
     imgui.Separator()
   imgui.Text("Files...")
 bottom_gui()


 end
end

graphite_main_window.add_module(my_dialog)
graphite_main_window.add_module(polycube_option_dialog)
































