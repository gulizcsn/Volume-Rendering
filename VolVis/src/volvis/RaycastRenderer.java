/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volvis;

import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2;
import com.jogamp.opengl.util.texture.Texture;
import com.jogamp.opengl.util.texture.awt.AWTTextureIO;

import gui.RaycastRendererPanel;
import gui.TransferFunction2DEditor;
import gui.TransferFunctionEditor;

import java.awt.image.BufferedImage;
import java.util.ArrayList;


import util.TFChangeListener;
import util.VectorMath;
import volume.GradientVolume;
import volume.Volume;
import volume.VoxelGradient;

/**
 * 
 * @author michel
 */
public class RaycastRenderer extends Renderer implements TFChangeListener {

	private Volume volume = null;
	private GradientVolume gradients = null;
	RaycastRendererPanel panel;
	TransferFunction tFunc;
	TransferFunctionEditor tfEditor;
	TransferFunction2DEditor tfEditor2D;

	private int rendererMode = 0;
	

	public RaycastRenderer() {
		panel = new RaycastRendererPanel(this);
		panel.setSpeedLabel("0");
	}

	public void setVolume(Volume vol) {
		System.out.println("Assigning volume");
		volume = vol;

		System.out.println("Computing gradients");
		gradients = new GradientVolume(vol);

		// set up image for storing the resulting rendering
		// the image width and height are equal to the length of the volume
		// diagonal
		int imageSize = (int) Math.floor(Math.sqrt(vol.getDimX()
				* vol.getDimX() + vol.getDimY() * vol.getDimY() + vol.getDimZ()
				* vol.getDimZ()));
		if (imageSize % 2 != 0) {
			imageSize = imageSize + 1;
		}
		image = new BufferedImage(imageSize, imageSize,
				BufferedImage.TYPE_INT_ARGB);
		// create a standard TF where lowest intensity maps to black, the
		// highest to white, and opacity increases
		// linearly from 0.0 to 1.0 over the intensity range
		tFunc = new TransferFunction(volume.getMinimum(), volume.getMaximum());

		// uncomment this to initialize the TF with good starting values for the
		// orange dataset
		tFunc.setTestFunc();

		tFunc.addTFChangeListener(this);
		tfEditor = new TransferFunctionEditor(tFunc, volume.getHistogram());

		tfEditor2D = new TransferFunction2DEditor(volume, gradients);
		tfEditor2D.addTFChangeListener(this);

		System.out.println("Finished initialization of RaycastRenderer");
	}

	public RaycastRendererPanel getPanel() {
		return panel;
	}

	public TransferFunction2DEditor getTF2DPanel() {
		return tfEditor2D;
	}

	public TransferFunctionEditor getTFPanel() {
		return tfEditor;
	}

        
        short getVoxel(double[] coord) {

            if (coord[0] < 0 || coord[0] > volume.getDimX() || coord[1] < 0 || coord[1] > volume.getDimY()
                    || coord[2] < 0 || coord[2] > volume.getDimZ()) {
                return 0;
            }

            int x = (int) Math.floor(coord[0]);
            int y = (int) Math.floor(coord[1]);
            int z = (int) Math.floor(coord[2]);

            return volume.getVoxel(x, y, z);
        }
        
        public short interpolationCalculation(double alfa, double beta, double gama, double c000,
			double c100, double c010, double c110, double c001, double c101,
			double c011, double c111){
            short result;

            result = (short) ((1-alfa)*(1-beta)*(1-gama)*c000+
                           (alfa)*(1-beta)*(1-gama)*c100+
                            (1-alfa)*(beta)*(1-gama)*c010+
                            (alfa)*(beta)*(1-gama)*c110+
                           (1-alfa)*(1-beta)*(gama)*c001+
                           (alfa)*(1-beta)*(gama)*c101+
                           (1-alfa)*(beta)*(gama)*c011+
                           (alfa)*(beta)*(gama)*c111);

             return result;
        
        }
        
	short voxelInterpolation(double[] coord) {

		if (coord[0] < 0 || coord[0] > volume.getDimX()-1 || coord[1] < 0
				|| coord[1] > volume.getDimY()-1 || coord[2] < 0
				|| coord[2] > volume.getDimZ()-1) {
			return 0;
		}
                
		int x = (int) Math.floor(coord[0]);
		int y = (int) Math.floor(coord[1]);
		int z = (int) Math.floor(coord[2]);
                
                if(x<0 || y<0 || z<0){
                
                    return 0;
                }
                
                int x1 = (int) Math.ceil(coord[0]);
		int y1 = (int) Math.ceil(coord[1]);
		int z1 = (int) Math.ceil(coord[2]);
              
                if(x1>volume.getDimX()-1 || y1>volume.getDimX()-1 || z1>volume.getDimX()-1){
                
                    return 0;
                } 
                
                if(this.interactiveMode)
                    {
                        return volume.getVoxel(x, y, z);
                    }        
        
                if(this.triLinear)
                    {
                        double alfa, beta, gama;
                        short interpolated;
                        alfa = coord[0] - x;
                        beta = coord[1] - y;
                        gama = coord[2] - z;

                        double c000 = volume.getVoxel(x, y, z);
                        double c100 = volume.getVoxel(x1, y, z);
                        double c010 = volume.getVoxel(x, y1, z);
                        double c110 = volume.getVoxel(x1, y1, z);
                        double c001 = volume.getVoxel(x, y, z1);
                        double c101 = volume.getVoxel(x1, y, z1);
                        double c011 = volume.getVoxel(x, y1, z1);
                        double c111 = volume.getVoxel(x1, y1, z1);


                        interpolated = interpolationCalculation(alfa, beta , gama , c000,
                                        c100, c010,c110,c001,c101,c011,c111);

                        return interpolated;
                    }
                
                else return volume.getVoxel(x, y, z);
              
	}

        short gradientInterpolation(double[] coord) {
                
		if (coord[0] < 0 || coord[0] > volume.getDimX()-1 || coord[1] < 0
				|| coord[1] > volume.getDimY()-1 || coord[2] < 0
				|| coord[2] > volume.getDimZ()-1) {
			return 0;
		}
                
		int x = (int) Math.floor(coord[0]);
		int y = (int) Math.floor(coord[1]);
		int z = (int) Math.floor(coord[2]);
                
                if(x<0 || y<0 || z<0){
                
                    return 0;
                }
                
                int x1 = (int) Math.ceil(coord[0]);
		int y1 = (int) Math.ceil(coord[1]);
		int z1 = (int) Math.ceil(coord[2]);
              
                if(x1>volume.getDimX()-1 || y1>volume.getDimX()-1 || z1>volume.getDimX()-1){
                
                    return 0;
                } 
                
                if(this.interactiveMode)
                    {
                        return (short)gradients.getGradient(x,y,z).mag;
                    }        
        
                      
                        double alfa, beta, gama;
                        short interpolatedGrad;
                        alfa = coord[0] - x;
                        beta = coord[1] - y;
                        gama = coord[2] - z;

                        double c000 = gradients.getGradient(x, y, z).mag;
                        double c100 = gradients.getGradient(x1, y, z).mag;
                        double c010 = gradients.getGradient(x, y1, z).mag;
                        double c110 = gradients.getGradient(x1, y1, z).mag;
                        double c001 = gradients.getGradient(x, y, z1).mag;
                        double c101 = gradients.getGradient(x1, y, z1).mag;
                        double c011 = gradients.getGradient(x, y1, z1).mag;
                        double c111 = gradients.getGradient(x1, y1, z1).mag;


                        interpolatedGrad = interpolationCalculation(alfa, beta , gama , c000,
                                        c100, c010,c110,c001,c101,c011,c111);
      
                        
                        return interpolatedGrad;
              
	}
        
        void setColor(int i, int j, int pixelColor,int stepSize)
        {
            if(stepSize==1)
            {
                image.setRGB(i, j, pixelColor);
            }
            else
            {
                for(int stepX = 0; stepX < stepSize; stepX++){
                    for(int stepY = 0; stepY < stepSize; stepY++){
                        if((j+stepX < image.getHeight()) &&(i+stepY < image.getWidth()))
                        {
                            image.setRGB(i+stepY,j+stepX,pixelColor);
                        }
                    }
                }
            }
        }
        
	void slicer(double[] viewMatrix) {

		// clear image
		for (int j = 0; j < image.getHeight(); j++) {
			for (int i = 0; i < image.getWidth(); i++) {
				image.setRGB(i, j, 0);
			}
		}
                   
		// vector uVec and vVec define a plane through the origin,
		// perpendicular to the view vector viewVec
		double[] viewVec = new double[3];
		double[] uVec = new double[3];
		double[] vVec = new double[3];
		VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6],
				viewMatrix[10]);
		VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
		VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

		// image is square
		int imageCenter = image.getWidth() / 2;

		double[] pixelCoord = new double[3];
		double[] volumeCenter = new double[3];
		VectorMath.setVector(volumeCenter, volume.getDimX() / 2,
				volume.getDimY() / 2, volume.getDimZ() / 2);

		// sample on a plane through the origin of the volume data
		double max = volume.getMaximum();
		TFColor voxelColor = new TFColor();

		for (int j = 0; j < image.getHeight(); j++) {
			for (int i = 0; i < image.getWidth(); i++) {
				pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0]
						* (j - imageCenter) + volumeCenter[0];
				pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1]
						* (j - imageCenter) + volumeCenter[1];
				pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2]
						* (j - imageCenter) + volumeCenter[2];

				int val = getVoxel(pixelCoord);
				

				// Map the intensity to a grey value by linear scaling
				voxelColor.r = val / max;
				voxelColor.g = voxelColor.r;
				voxelColor.b = voxelColor.r;
				voxelColor.a = val > 0 ? 1.0 : 0.0; // this makes intensity 0
													// completely transparent
													// and the rest opaque
				// Alternatively, apply the transfer function to obtain a color
				// voxelColor = tFunc.getColor(val);

				// BufferedImage expects a pixel color packed as ARGB in an int
				int c_alpha = voxelColor.a <= 1.0 ? (int) Math
						.floor(voxelColor.a * 255) : 255;
				int c_red = voxelColor.r <= 1.0 ? (int) Math
						.floor(voxelColor.r * 255) : 255;
				int c_green = voxelColor.g <= 1.0 ? (int) Math
						.floor(voxelColor.g * 255) : 255;
				int c_blue = voxelColor.b <= 1.0 ? (int) Math
						.floor(voxelColor.b * 255) : 255;
				int pixelColor = (c_alpha << 24) | (c_red << 16)
						| (c_green << 8) | c_blue;
				image.setRGB(i, j, pixelColor);
			}
		}

	}

	void mip(double[] viewMatrix) {

		// clear image
		for (int j = 0; j < image.getHeight(); j++) {
			for (int i = 0; i < image.getWidth(); i++) {
				image.setRGB(i, j, 0);
			}
		}

		// vector uVec and vVec define a plane through the origin,
		// perpendicular to the view vector viewVec
		double[] viewVec = new double[3];
		double[] uVec = new double[3];
		double[] vVec = new double[3];
		VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6],viewMatrix[10]);
		VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
		VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

		// image is square
		int imageCenter = image.getWidth() / 2;

		double[] pixelCoord = new double[3];
		
		double[] volumeCenter = new double[3];
		VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

		TFColor voxelColor = new TFColor();

		int[] volumeDimensions = new int[] { volume.getDimX(), volume.getDimY(), volume.getDimZ() };
		
		// using the maximum volume dimension to bound the ray. Can be optimized. 
		int maxDimension = getMax(volumeDimensions);
		int minDimension = -1 * maxDimension; 
		
		// When interactive mode, number of samples taken is reduced to make program responsive and fast.
		int stepSize = 1;
		if(this.interactiveMode){
			stepSize = 3;
		}
		double max = volume.getMaximum();
		
		int maximum = 0;
		for (int j = 0; j < image.getHeight(); j+=stepSize) {
			for (int i = 0; i < image.getWidth(); i+=stepSize) {

				int maxVal = 0;
				for( double k=minDimension; k<= maxDimension; k++){
					
					pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] 
						* (j - imageCenter) + volumeCenter[0] + k * viewVec[0];
					pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1]
						* (j - imageCenter) + volumeCenter[1] + k
						* viewVec[1];
					pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2]
						* (j - imageCenter) + volumeCenter[2] + k
						* viewVec[2];
				
				
				int val = voxelInterpolation(pixelCoord);
				
				
				if (val > maxVal)
					maxVal = val;

				}
				
			    if(maxVal > maximum){
			    	maximum = maxVal;
			    }
				// Map the intensity to a grey value by linear scaling
				 voxelColor.r = maxVal/max;
				 voxelColor.g = voxelColor.r;
				 voxelColor.b = voxelColor.r;
				 voxelColor.a = maxVal > 0 ? 1.0 : 0.0; // this makes
				// intensity 0 completely transparent and the rest opaque
				// Alternatively, apply the transfer function to obtain a
				// color
				//voxelColor = tFunc.getColor(maxVal);
				// BufferedImage expects a pixel color packed as ARGB in an int
				int c_alpha = voxelColor.a <= 1.0 ? (int) Math
						.floor(voxelColor.a * 255) : 255;
				int c_red = voxelColor.r <= 1.0 ? (int) Math
						.floor(voxelColor.r * 255) : 255;
				int c_green = voxelColor.g <= 1.0 ? (int) Math
						.floor(voxelColor.g * 255) : 255;
				int c_blue = voxelColor.b <= 1.0 ? (int) Math
						.floor(voxelColor.b * 255) : 255;
				int pixelColor = (c_alpha << 24) | (c_red << 16)
						| (c_green << 8) | c_blue;
				setColor(i,j,pixelColor,stepSize);
                               // image.setRGB(i, j, pixelColor);
			}
		}
		
		
	}

	
	 void compositing(double[] viewMatrix, boolean shade) {

	        // clear image
	        for (int j = 0; j < image.getHeight(); j++) {
	            for (int i = 0; i < image.getWidth(); i++) {
	                image.setRGB(i, j, 0);
	            }
	        }
	       
	        // vector uVec and vVec define a plane through the origin, 
	        // perpendicular to the view vector viewVec
	        double[] viewVec = new double[3];
	        double[] uVec = new double[3];
	        double[] vVec = new double[3];
	        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
	        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
	        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);
	        
	        // image is square
	        int imageCenter = image.getWidth() / 2;

	        double[] pixelCoord = new double[3];
	        double[] volumeCenter = new double[3];
	        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);
	        
	        int[] volumeDimensions = new int[] {volume.getDimX(), volume.getDimY(), volume.getDimZ()};
	        int maxDimension = getMax(volumeDimensions);
	        
	        // sample on a plane through the origin of the volume data
	        double max = volume.getMaximum();
	        ArrayList<TFColor> colors = new ArrayList<TFColor>(3);
	        int stepSize=1;
            if(this.interactiveMode){
            	//Working with lower resolution to make the application's response time better.
            	stepSize = 3;
            }
	        for (int j = 0; j < image.getHeight(); j+=stepSize) {
	            for (int i = 0; i < image.getWidth(); i+=stepSize) {
	                
	                colors.add(0, new TFColor());
	                colors.add(1,new TFColor());
	                colors.add(2, new TFColor());
	                
	                
	                for (double t=(-1*maxDimension); t<= maxDimension; t++) {
	                // Optimization possible by step size
	                    
	                    pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
	                            + volumeCenter[0] + t * viewVec[0];
	                    pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
	                            + volumeCenter[1] + t * viewVec[1];
	                    pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
	                            + volumeCenter[2] + t * viewVec[2];
                            
                            int val = voxelInterpolation(pixelCoord);
                            TFColor voxelColor = tFunc.getColor(val);
                            
                            if(shade)
                                voxelColor = phongShading(voxelColor,viewVec,0.1,0.7,0.2,10,pixelCoord);
                                
	                    colors = getDVRColor(voxelColor,colors.get(1),colors.get(2));
	                }
	                
	                int c_alpha = (1 - colors.get(1).a) <= 1.0 ? (int) Math.floor((1 - colors.get(1).a) * 255) : 255;
	                int c_red = colors.get(1).r <= 1.0 ? (int) Math.floor(colors.get(1).r * 255) : 255;
	                int c_green = colors.get(1).g <= 1.0 ? (int) Math.floor(colors.get(1).g * 255) : 255;
	                int c_blue = colors.get(1).b <= 1.0 ? (int) Math.floor(colors.get(1).b * 255) : 255;
	                    
	                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
	                setColor(i,j,pixelColor,stepSize);
                        //image.setRGB(i, j, pixelColor);
	            }
	        } 
	    }
	   
        void opacityWeighting(double[] viewMatrix, boolean shade) {

// clear image
		for (int j = 0; j < image.getHeight(); j++) {
			for (int i = 0; i < image.getWidth(); i++) {
				image.setRGB(i, j, 0);
			}
		}
                short intensity = tfEditor2D.triangleWidget.baseIntensity;
                double radius = tfEditor2D.triangleWidget.radius;
                TFColor color = tfEditor2D.triangleWidget.color;
                double maxGradientThreshold = tfEditor2D.triangleWidget.maxGradientThreshold;
                double minGradientThreshold = tfEditor2D.triangleWidget.minGradientThreshold;
               
                
                short interpolatedgradmag;
		// vector uVec and vVec define a plane through the origin,
		// perpendicular to the view vector viewVec
		double[] viewVec = new double[3];
		double[] uVec = new double[3];
		double[] vVec = new double[3];
		VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6],viewMatrix[10]);
		VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
		VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

		// image is square
		int imageCenter = image.getWidth() / 2;

		double[] pixelCoord = new double[3];
		
		double[] volumeCenter = new double[3];
		VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

		TFColor voxelColor = new TFColor();

		int[] volumeDimensions = new int[] { volume.getDimX(), volume.getDimY(), volume.getDimZ() };
		
		// using the maximum volume dimension to bound the ray. Can be optimized. 
		int maxDimension = getMax(volumeDimensions);
		int minDimension = -1 * maxDimension; 
		
		// When interactive mode, number of samples taken is reduced to make program responsive and fast.
		int stepSize = 1;
                ArrayList<TFColor> colors = new ArrayList<TFColor>(3);
		if(this.interactiveMode){
			stepSize = 3;
		}
		for (int j = 0; j < image.getHeight(); j+=stepSize) {
			for (int i = 0; i < image.getWidth(); i+=stepSize) {
                            colors.add(0, new TFColor());
                            colors.add(1,new TFColor());
                            colors.add(2, new TFColor());
	                
			
				for( double k=minDimension; k<= maxDimension; k++){
					
					pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] 
						* (j - imageCenter) + volumeCenter[0] + k * viewVec[0];
					pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1]
						* (j - imageCenter) + volumeCenter[1] + k
						* viewVec[1];
					pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2]
						* (j - imageCenter) + volumeCenter[2] + k
						* viewVec[2];
                                        
                                        
                                        int val = voxelInterpolation(pixelCoord);
                                        voxelColor.r = color.r;
                                        voxelColor.g = color.g;
                                        voxelColor.b = color.b;
                              
                               
                                        interpolatedgradmag = gradientInterpolation(pixelCoord);
                                    
                                        if(interpolatedgradmag==0 && val==intensity)
                                        {
                                           voxelColor.a = color.a;
                                           
                                           if(shade)
                                                voxelColor = phongShading(voxelColor,viewVec,0.1,0.7,0.2,10,pixelCoord);
    
                                        }
                                        else if((interpolatedgradmag > 0) && (interpolatedgradmag < maxGradientThreshold) && (interpolatedgradmag > minGradientThreshold) && ((val-radius*interpolatedgradmag) <= intensity) && ((val+radius*interpolatedgradmag) >= intensity) )
                                        {
                                           
                                           voxelColor.a = color.a * (1.0 - (1/radius)*(double)Math.abs((intensity-val)/interpolatedgradmag)); 
                                           
                                           if(shade)
                                                voxelColor = phongShading(voxelColor,viewVec,0.1,0.7,0.2,10,pixelCoord);      
                                        }
                                        else {
                                           voxelColor.a = 0;
                                        }
                                 
                                        colors = getDVRColor(voxelColor,colors.get(1),colors.get(2));
                               }
                                
                                TFColor shadedColor = colors.get(1);
				// BufferedImage expects a pixel color packed as ARGB in an int
				int c_alpha = (1 - shadedColor.a) <= 1.0 ? (int) Math.floor((1 - shadedColor.a) * 255) : 255;
                                int c_red = shadedColor.r <= 1.0 ? (int) Math.floor(shadedColor.r * 255) : 255;
                                int c_green = shadedColor.g <= 1.0 ? (int) Math.floor(shadedColor.g * 255) : 255;
                                int c_blue = shadedColor.b <= 1.0 ? (int) Math.floor(shadedColor.b * 255) : 255;
                                
                                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                                setColor(i,j,pixelColor,stepSize);
                                //image.setRGB(i, j, pixelColor);
              
			}
		}
    }

        public ArrayList<TFColor> getDVRColor(TFColor voxelColor, TFColor nextColor, TFColor prevColor) {
            ArrayList<TFColor> colors = new ArrayList<TFColor>(3);
            
            nextColor.r = voxelColor.a * voxelColor.r + (1 - voxelColor.a) * prevColor.r;
            nextColor.g = voxelColor.a * voxelColor.g + (1 - voxelColor.a) * prevColor.g;
            nextColor.b = voxelColor.a * voxelColor.b + (1 - voxelColor.a) * prevColor.b;
                
            nextColor.a = (1 - voxelColor.a) * prevColor.a;

            prevColor = nextColor;

            //Return voxelColor, nextColor and prevColor 
            colors.add(voxelColor);
            colors.add(nextColor);
            colors.add(prevColor);

            return colors;
        
	}
 
        public TFColor phongShading(TFColor color,double [] viewVec, double kambient, double kdiff, double kspec, double shine, double pixelCoord[])
        {
            TFColor newColor = new TFColor();
            if (pixelCoord[0] < 0 || pixelCoord[0] > volume.getDimX() - 1|| pixelCoord[1] < 0 || pixelCoord[1] > volume.getDimY() - 1
                || pixelCoord[2] < 0 || pixelCoord[2] > volume.getDimZ() - 1) {
                return color;
            }
            
            VoxelGradient voxelgrad = gradients.getGradient((int)Math.floor(pixelCoord[0]),(int)Math.floor(pixelCoord[1]),(int)Math.floor(pixelCoord[2]));
           
            //normalized view Vector
            double normalview [] = new double [] {-viewVec[0]/VectorMath.length(viewVec),-viewVec[1]/VectorMath.length(viewVec),-viewVec[2]/VectorMath.length(viewVec)};
            //normalized grad vector
            double grad[] = new double [] {-voxelgrad.x/voxelgrad.mag,-voxelgrad.y/voxelgrad.mag,-voxelgrad.z/voxelgrad.mag};
            
            //diffuse weighting
            //grad vector is used as an estimation of N (normal vector)
            double diffuseproduct = VectorMath.dotproduct(normalview, grad);
            TFColor color2 = new TFColor(color.r, color.g, color.b, color.a);
            if(diffuseproduct>0){
                
                color2.r = kambient + color.r * kdiff * diffuseproduct;
                color2.g = kambient + color.g * kdiff * diffuseproduct;
                color2.b = kambient + color.b * kdiff * diffuseproduct;
            }
            
            //spec weighting
            //grad vector is used as an estimation of N (normal vector)
            //as we assume L=V, halfway vector is V, which is view vector
            
            double specCal= VectorMath.dotproduct(grad,normalview);
             if(specCal>0){
                
                double specproduct = Math.pow(specCal, shine);

                color2.r += kspec * specproduct;
                color2.g += kspec * specproduct;
                color2.b += kspec * specproduct;
            }

            newColor.a = color2.a;
            newColor.r = color2.r;
            newColor.g = color2.g;
            newColor.b = color2.b;
            return newColor;
        }

	private void drawBoundingBox(GL2 gl) {
		gl.glPushAttrib(GL2.GL_CURRENT_BIT);
		gl.glDisable(GL2.GL_LIGHTING);
		gl.glColor4d(1.0, 1.0, 1.0, 1.0);
		gl.glLineWidth(1.5f);
		gl.glEnable(GL.GL_LINE_SMOOTH);
		gl.glHint(GL.GL_LINE_SMOOTH_HINT, GL.GL_NICEST);
		gl.glEnable(GL.GL_BLEND);
		gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

		gl.glBegin(GL.GL_LINE_LOOP);
		gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0,
				volume.getDimZ() / 2.0);
		gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0,
				volume.getDimZ() / 2.0);
		gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0,
				volume.getDimZ() / 2.0);
		gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0,
				volume.getDimZ() / 2.0);
		gl.glEnd();

		gl.glBegin(GL.GL_LINE_LOOP);
		gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0,
				-volume.getDimZ() / 2.0);
		gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0,
				-volume.getDimZ() / 2.0);
		gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0,
				-volume.getDimZ() / 2.0);
		gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0,
				-volume.getDimZ() / 2.0);
		gl.glEnd();

		gl.glBegin(GL.GL_LINE_LOOP);
		gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0,
				-volume.getDimZ() / 2.0);
		gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0,
				volume.getDimZ() / 2.0);
		gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0,
				volume.getDimZ() / 2.0);
		gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0,
				-volume.getDimZ() / 2.0);
		gl.glEnd();

		gl.glBegin(GL.GL_LINE_LOOP);
		gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0,
				-volume.getDimZ() / 2.0);
		gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0,
				volume.getDimZ() / 2.0);
		gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0,
				volume.getDimZ() / 2.0);
		gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0,
				-volume.getDimZ() / 2.0);
		gl.glEnd();

		gl.glBegin(GL.GL_LINE_LOOP);
		gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0,
				-volume.getDimZ() / 2.0);
		gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0,
				volume.getDimZ() / 2.0);
		gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0,
				volume.getDimZ() / 2.0);
		gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0,
				-volume.getDimZ() / 2.0);
		gl.glEnd();

		gl.glBegin(GL.GL_LINE_LOOP);
		gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0,
				-volume.getDimZ() / 2.0);
		gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0,
				volume.getDimZ() / 2.0);
		gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0,
				volume.getDimZ() / 2.0);
		gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0,
				-volume.getDimZ() / 2.0);
		gl.glEnd();

		gl.glDisable(GL.GL_LINE_SMOOTH);
		gl.glDisable(GL.GL_BLEND);
		gl.glEnable(GL2.GL_LIGHTING);
		gl.glPopAttrib();

	}

	// Method for getting the maximum value
	public static int getMax(int[] inputArray) {

		int maxValue = inputArray[0];

		for (int i = 1; i < inputArray.length; i++) {

			if (inputArray[i] > maxValue) {

				maxValue = inputArray[i];
			}
		}

		return maxValue;
	}

	@Override
	public void visualize(GL2 gl) {

		if (volume == null) {
			return;
		}

		drawBoundingBox(gl);

		gl.glGetDoublev(GL2.GL_MODELVIEW_MATRIX, viewMatrix, 0);

		long startTime = System.currentTimeMillis();

		// rendererModes: 0 = slicer, 1 = mip, 2 = compositing, 3 = 2dtransfer..
		switch (rendererMode) {
                case 0:
			slicer(viewMatrix);
			break;
		case 1:
			mip(viewMatrix);
			break;
		case 2:
			compositing(viewMatrix,false);
			break;
                case 3:
			opacityWeighting(viewMatrix,false);
			break;
                case 4:
			compositing(viewMatrix,true);
			break;
                case 5:
			opacityWeighting(viewMatrix,true);
			break;        
		default:
			System.out.println("Invalid renderer method.");
			break;
		}

		long endTime = System.currentTimeMillis();
		double runningTime = (endTime - startTime);
		panel.setSpeedLabel(Double.toString(runningTime));

		Texture texture = AWTTextureIO.newTexture(gl.getGLProfile(), image,
				false);

		gl.glPushAttrib(GL2.GL_LIGHTING_BIT);
		gl.glDisable(GL2.GL_LIGHTING);
		gl.glEnable(GL.GL_BLEND);
		gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

		// draw rendered image as a billboard texture
		texture.enable(gl);
		texture.bind(gl);
		double halfWidth = image.getWidth() / 2.0;
		gl.glPushMatrix();
		gl.glLoadIdentity();
		gl.glBegin(GL2.GL_QUADS);
		gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
		gl.glTexCoord2d(0.0, 0.0);
		gl.glVertex3d(-halfWidth, -halfWidth, 0.0);
		gl.glTexCoord2d(0.0, 1.0);
		gl.glVertex3d(-halfWidth, halfWidth, 0.0);
		gl.glTexCoord2d(1.0, 1.0);
		gl.glVertex3d(halfWidth, halfWidth, 0.0);
		gl.glTexCoord2d(1.0, 0.0);
		gl.glVertex3d(halfWidth, -halfWidth, 0.0);
		gl.glEnd();
		texture.disable(gl);
		texture.destroy(gl);
		gl.glPopMatrix();

		gl.glPopAttrib();

		if (gl.glGetError() > 0) {
			System.out.println("some OpenGL error: " + gl.glGetError());
		}

	}

	private BufferedImage image;
	private double[] viewMatrix = new double[4 * 4];
        private boolean triLinear;
	@Override
	public void changed() {
		for (int i = 0; i < listeners.size(); i++) {
			listeners.get(i).changed();
		}
	}

	public void setMode(int mode) {
		rendererMode = mode;
                changed();
	}
        
        public int getMode() {
		return rendererMode;
	}
        
        public void setTriLinear(boolean triLinear) {
			this.triLinear = triLinear;
			 changed();
	}
}
