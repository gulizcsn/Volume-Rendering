/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volume;

/**
 *
 * @author michel
 */
public class GradientVolume {

    public GradientVolume(Volume vol) {
        volume = vol;
        dimX = vol.getDimX();
        dimY = vol.getDimY();
        dimZ = vol.getDimZ();
        data = new VoxelGradient[dimX * dimY * dimZ];
        compute();
        maxmag = -1.0;
    }

    public VoxelGradient getGradient(int x, int y, int z) {
        return data[x + dimX * (y + dimY * z)];
    }

    
    public void setGradient(int x, int y, int z, VoxelGradient value) {
        data[x + dimX * (y + dimY * z)] = value;
    }

    public void setVoxel(int i, VoxelGradient value) {
        data[i] = value;
    }

    public VoxelGradient getVoxel(int i) {
        return data[i];
    }

    public int getDimX() {
        return dimX;
    }

    public int getDimY() {
        return dimY;
    }

    public int getDimZ() {
        return dimZ;
    }

 private void compute() {

        // this just initializes all gradients to the vector (0,0,0)
        System.out.println(data.length);
        for (int i=0; i<data.length; i++) {
            data[i] = zero;
        }
        
        int dx=1;
        int dy=1;
        int dz=1;
        float voxelgrad_x,voxelgrad_y,voxelgrad_z;
<<<<<<< Updated upstream:VolVis/src/volume/GradientVolume.java
       // to take get the values between voxel range we need to calculate dimensions exceed to this range or not
//       int rangeforX=Math.min(getDimX(),volume.getMaximum());
//       int rangeforY=Math.min(getDimY(),volume.getMaximum());
//       int rangeforZ=Math.min(getDimZ(),volume.getMaximum());
=======

>>>>>>> Stashed changes:visualization-master-3/VolVis/src/volume/GradientVolume.java
        int rangeforX=volume.getDimX();
        int rangeforY=volume.getDimY();
        int rangeforZ=volume.getDimZ();
     
        for (int j=0; j<rangeforX; j++) {
            for (int k=0; k<rangeforY; k++) {
                for (int l=0; l<rangeforZ; l++) {
                    
                // gradient of the outer frame of pixels should be zero
                     if((j==0)||(k==0)||(l==0)||(j==(rangeforX-1))||(k==(rangeforY-1))||(l==(rangeforZ-1)))
                     {
                        this.setGradient(j, k, l, zero);
                     }
                     else
                     {
               //derivative approximation for gradient
                        voxelgrad_x= (volume.getVoxel(j+dx, k, l) - volume.getVoxel(j-dx, k, l))/2.0f;
                        voxelgrad_y= (volume.getVoxel(j, k+dy, l) - volume.getVoxel(j, k-dy, l))/2.0f;
                        voxelgrad_z= (volume.getVoxel(j, k, l+dz) - volume.getVoxel(j, k, l-dz))/2.0f;
                        VoxelGradient gradient = new VoxelGradient(voxelgrad_x,voxelgrad_y,voxelgrad_z); 
                        this.setGradient(j,k,l,gradient);
                     }
                }
            }
        }
    }

    
    public double getMaxGradientMagnitude() {
        if (maxmag >= 0) {
            return maxmag;
        } else {
            double magnitude = data[0].mag;
            for (int i=0; i<data.length; i++) {
                magnitude = data[i].mag > magnitude ? data[i].mag : magnitude;
            }   
            maxmag = magnitude;
            return magnitude;
        }
    }
    
    private int dimX, dimY, dimZ;
    private VoxelGradient zero = new VoxelGradient();
    VoxelGradient[] data;
    Volume volume;
    double maxmag;
}
