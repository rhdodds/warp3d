:: Delete obj and mod files
for %%F in (8to20 add_elm mesh_plot mesh_ssy mesh_ssy2_ts pipe_mesh_gen mesh_scp mesh_cell) do (
cd %%F
del *.obj /q 
del *.mod /q 
cd ..
)