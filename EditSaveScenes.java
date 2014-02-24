// STAR-CCM+ macro
// This runs a STAR-CCM+ job in parallel under a batch queing system
// R. Reynolds
//=======================================
package macro;

import java.util.*;
import star.flow.*;
import star.base.neo.*;   // for use with version >= 1.09
import star.base.report.*;
import star.common.*;
import star.coremodule.services.*;
import star.saturb.*;
import star.keturb.*;
import star.kwturb.*;
//import star.scenefile.*;
import star.common.StarPlot.*;
import star.vis.*;
//import star.segregatedenergy.*;
//import star.segregatedflow.*;
import star.coupledflow.*;
import star.meshing.*;



public class EditSaveScenes extends StarMacro {

  public void execute() {

    // Get the simulation in preparation for running
    Simulation sim = getActiveSimulation();
    String currentDirectory = sim.getSessionDir();

    //Setup Alpha and Beta

    Region region_0 =
      sim.getRegionManager().getRegion("Domain");

    // Adjusts the Wind Axis Coordinate System for the AOA
    //Get Volume Mesh:
    FvRepresentation volumeMesh = 
      ((FvRepresentation) sim.getRepresentationManager().getObject("Volume Mesh"));
    //Get vector of all scenes:
    Collection<Scene> colSCN = sim.getSceneManager().getScenes();
        
    if (!colSCN.isEmpty()) {//Make sure scenes exist
        for (Scene sce : colSCN){//Save all scenes
            CurrentView currentView = 
            sce.getCurrentView();
            sce.getDisplayerManager().setRepresentation(volumeMesh);
            sce.export3DSceneFileAndWait(resolvePath(currentDirectory+"\\"+sce+".sce"), true);
        }
    }
    Report r = sim.getReportManager().getReport(s);
    //r.getReferences().toString();
    r.printReport(null, true);
    Collection<MonitorReport> colMR= sim.getReportManager().getObjectsOf(Collection<ReportMonitor>)
            for (MonitorReport mr : colMR) {
                sim.println(mr.);
            }
    //String PlotNameVector[] ={"CDBA Monitor Plot", "CLBA Monitor Plot"};
    //for (String s2 : PlotNameVector){
      //  MonitorPlot plot = ((MonitorPlot) sim.getPlotManager().getObject(s2));
        //plot.encode(resolvePath(currentDirectory+"\\"+s2+".png"), "png", 1920, 1080);
    //}
        
  }
}





    
