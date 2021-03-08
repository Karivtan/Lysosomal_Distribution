package LeidenUniv.Fish;
import ij.*;
import ij.process.*;
import ij.gui.*;
import ij.plugin.*;
import ij.plugin.frame.*;
import ij.measure.*;

public class Lysosomal_distribution implements PlugIn {

	public void run(String arg) {
		new WaitForUserDialog("Select the lysosome image").show();
		ImagePlus lyso = IJ.getImage();
		
		new WaitForUserDialog("Select the Macrophage image").show();
		ImagePlus macro = IJ.getImage();
		Calibration cal = macro.getCalibration();
		
		RoiManager rm = RoiManager.getInstance();
		if (rm!=null){
			rm.close();
			rm = new RoiManager();
		} else {
			rm = new RoiManager();
		}
		ImagePlus imp2 = macro.duplicate();
		IJ.run(imp2, "Gaussian Blur...", "sigma=1 stack");
		IJ.run(imp2, "Convert to Mask", "method=Li background=Dark calculate black");
		IJ.run(imp2, "Analyze Particles...", "size=15-Infinity show=Nothing add stack");

		Roi [] imageRois = rm.getRoisAsArray();
		IJ.setForegroundColor(0, 0, 0);
		//imp2.show();
		int nRois=rm.getCount();
		ResultsTable rt = new ResultsTable();
		for (int i=0; i<nRois;i++){
			rt.incrementCounter();
			Roi cRoi =imageRois[i]; // curent roi
			macro.setRoi(cRoi);
			imp2.setSlice(cRoi.getZPosition());
			ImageStatistics is = macro.getStatistics(Measurements.CENTROID+Measurements.ELLIPSE); // get statistics
			
			double cAngle= is.angle;
			double radAngle =(cAngle+90)/180*Math.PI;
			double mradAngle =(cAngle/180)*Math.PI;
			double cCos= Math.cos(radAngle);
			double cSin= Math.sin(radAngle);
			double ccCos= Math.cos(mradAngle);
			double ccSin= Math.sin(mradAngle);
			double unit = cal.pixelHeight;
			double cx=cCos*((is.minor/unit)+1);
			double cy=cSin*((is.minor/unit)+1);
			double ccx=ccCos*((is.major/unit)+1);
			double ccy=ccSin*((is.major/unit)+1);
			int[] xpoints = {(int)(is.xCentroid/unit-cx),(int)(is.xCentroid/unit+cx),(int)(is.xCentroid/unit+cx+ccx),(int)(is.xCentroid/unit-cx+ccx)};
			int[] ypoints = {(int)(is.yCentroid/unit+cy),(int)(is.yCentroid/unit-cy),(int)(is.yCentroid/unit-cy-ccy),(int)(is.yCentroid/unit+cy-ccy)};
			imp2.setRoi(new PolygonRoi(xpoints,ypoints,4,Roi.POLYGON));
			rm.addRoi(imp2.getRoi());
			rm.setSelectedIndexes(new int[]{i,(2*i+nRois)});
			rm.runCommand(imp2,"AND");
			Roi topRoi = imp2.getRoi();
			macro.setSlice(cRoi.getZPosition());
			macro.setRoi(topRoi);
			lyso.setSlice(cRoi.getZPosition());
			lyso.setRoi(topRoi);
			double topM = macro.getStatistics(2).mean;
			double topL = lyso.getStatistics(2).mean;
			
			int[] xpoints2 = {(int)(is.xCentroid/unit-cx),(int)(is.xCentroid/unit+cx),(int)(is.xCentroid/unit+cx-ccx),(int)(is.xCentroid/unit-cx-ccx)};
			int[] ypoints2 = {(int)(is.yCentroid/unit+cy),(int)(is.yCentroid/unit-cy),(int)(is.yCentroid/unit-cy+ccy),(int)(is.yCentroid/unit+cy+ccy)};
			imp2.setRoi(new PolygonRoi(xpoints2,ypoints2,4,Roi.POLYGON));
			rm.addRoi(imp2.getRoi());
			rm.setSelectedIndexes(new int[]{i,(2*i+nRois+1)});
			rm.runCommand(imp2,"AND");
			Roi bottomRoi = imp2.getRoi();
			macro.setRoi(bottomRoi);
			lyso.setRoi(bottomRoi);
			double botM = macro.getStatistics(2).mean;
			double botL = lyso.getStatistics(2).mean;
			rt.addValue("X position", is.xCentroid);
			rt.addValue("Y position", is.yCentroid);
			rt.addValue("Z position", cRoi.getZPosition());
			
			rt.addValue("Top Macrophage", topM);
			rt.addValue("Bottom Macrophage", botM);
			rt.addValue("Ratio Top/Bottom Macrophages", topM/botM);
			rt.addValue("Top Lysozomes", topL);
			rt.addValue("Bottom Lysozomes", botL);
			rt.addValue("Ratio Top/Bottom Lysozomes", topL/botL);
			rt.addValue("Percentage Top", topL/(topL+botL));
			rt.addValue("Percentage Bottom", botL/(topL+botL));
			
		}
		rt.show("Mean intensities of top and bottom of fitted ellipses");
	}
}