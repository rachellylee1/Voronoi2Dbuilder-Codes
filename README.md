# Voronoi2Dbuilder-Codes
 2DRealisticMicrostructure Builder

1. Make sure you have the following files in the same directory.  “Voronoi2Dbuilder”, “cross2D”, and “VoronoiDXFApplication.swp” 
 
2. Open “Voronoi2Dbuilder” 
 
3. Scroll down to line 358 of the code:  fileID = fopen(sprintf('G:\\Vertex\\Vertex%d.txt',icell),'w'); 
 
4. Replace, G:, with the directory of your choice. 
 
5. Create a folder in that directory and name it “Vertex”
 
6. Into the Command Window, type in Voronoi2Dbuilder(ngrains, bndwidth, mult, bins) 
 
7. Replace ngrains, bndwidth, mult, and bins with the number of grains, grain boundary width, average grain diameter, and number of bins for the histogram respectively. Remember to assume the program is already in um scale. Ex. I want to build a 50 grain sample with 1 nm grain boundary width. The average grain size is 100 µm. I want the histogram of grain sizes to have 5 bins. ngrains = 50, bndwidth = 0.001, mult = 100, and bins = 5 Type in to the Command Window: Voronoi2Dbuilder(50, 0.001, 100, 5) 
 
8. Open Solidworks. 
 
9. Make a new part (at top of the screen)
Make sure “Part” is highlighted and press “OK” 
 
10. Go to “Tools” at the top and scroll down to “ Options…” and select. 
 
11. Go to “Relations/Snaps” on the left hand side. 
 
12. Disable “Enable snapping” and press “OK”. 
 
13. Go back to “Tools” and select “Macro" “Edit…” 
 
14. Find the “VoronoiDXFApplication.swp” file and open. 
 
15. In the pop up window, find and press the small green arrow button towards the top of the window.  The realistic model should automatically build in Solidworks. 
 
16. Name the file and save as a .dxf 
 
17. On the left hand side of the Solidworks window, make sure “Annotation view” and “*Front” is selected. Unselect other views from “Views to Export”. 
 
18. Press the green check mark. 
 
19. Save the DXF file. 
 
 
