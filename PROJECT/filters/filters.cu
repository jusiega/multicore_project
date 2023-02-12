#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 199309L
#endif
#include <sys/time.h>
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <unistd.h>


#define STB_IMAGE_IMPLEMENTATION
#include "../Common/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../Common/stb_image_write.h"








double wtime()
{
static int sec = -1;
struct timespec tv;
clock_gettime(CLOCK_MONOTONIC, &tv);
if (sec < 0) sec = tv.tv_sec;
return (tv.tv_sec - sec) + 1.0e-9*tv.tv_nsec;
}



__global__ void fblur(float* r, float* g, float* b, float* r2, float* g2, float* b2, int width, int height, float*  blur, int blursize) //standard blur
{   
    __syncthreads();
    long long int idx= threadIdx.x + blockIdx.x * blockDim.x;
    int j=int(idx%width);
    int i=int(idx/width);

    float red_tmp=0;
    float green_tmp=0;
    float blu_tmp=0;

    int shift=(blursize-1)/2;
__syncthreads();
    //MAIN BODY
    if(i>=shift && i<height-shift && j>=shift && j<width-shift)
    {

            for (int k=0; k<blursize; k++)
            {
                for (int l=0; l<blursize; l++)
                {
                    red_tmp+=blur[k*blursize+l]*(r[(i-shift+k)*width + (j-shift+l)]);
                    green_tmp+=blur[k*blursize+l]*(g[(i-shift+k)*width + (j-shift+l)]);
                    blu_tmp+=blur[k*blursize+l]*(b[(i-shift+k)*width + (j-shift+l)]);
                    //blu_tmp+=blur[k*blursize+l]*float(b[i-(blursize-1)/2+k][j-(blursize-1)/2+l]);
                }
                
            }
            r2[i*width+j]=(red_tmp);
            g2[i*width+j]=(green_tmp);
            b2[i*width+j]=(blu_tmp);
            red_tmp=0;
            green_tmp=0;
            blu_tmp=0;
    }
    //EDGE HANDLING HERE IS BROKEN SO IT HAS BEEN TURNED OFF
    //EDG TOP
    /*__syncthreads();
    if(i>=0 && i<shift && j>=0 && j<width)
    {
            for (int k=0; k<blursize; k++)
            {
                for (int l=0; l<blursize; l++)
                {
                    red_tmp+=blur[k*blursize+l]*r[abs(((i-shift)+k))*width + abs((j-shift+l))];
                    green_tmp+=blur[k*blursize+l]*g[abs(((i-shift)+k))*width + abs((j-shift+l))];
                    blu_tmp+=blur[k*blursize+l]*b[abs(((i-shift)+k))*width + abs((j-shift+l))];
                    
                }
            }
            r2[i*width+j]=(red_tmp);
            g2[i*width+j]=(green_tmp);
            b2[i*width+j]=(blu_tmp);
            red_tmp=0;
            green_tmp=0;
            blu_tmp=0;
    }
    __syncthreads();
    //EDG left
    if(i>=0 && i<height && j>=0 && j<shift)
    {
            for (int k=0; k<blursize; k++)
            {
                for (int l=0; l<blursize; l++)
                {
                    red_tmp+=blur[k*blursize+l]*r[abs(((i-shift)+k))*width + abs((j-shift+l))];
                    green_tmp+=blur[k*blursize+l]*g[abs(((i-shift)+k))*width + abs((j-shift+l))];
                    blu_tmp+=blur[k*blursize+l]*b[abs(((i-shift)+k))*width + abs((j-shift+l))];
                    
                }
            }
            r2[i*width+j]=(red_tmp);
            g2[i*width+j]=(green_tmp);
            b2[i*width+j]=(blu_tmp);
            red_tmp=0;
            green_tmp=0;
            blu_tmp=0;
    }
    __syncthreads();
    //EDGE DOWN
    if(i>=height-shift && i<height && j>=0 && j<width)
    {
            for (int k=0; k<blursize; k++)
            {
                for (int l=0; l<blursize; l++)
                {
                    red_tmp+=blur[k*blursize+l]*(r[(-1+height-abs(height-1-(i-shift+k)))*width + width-1-abs(width-1-(j-shift+l))]);
                    green_tmp+=blur[k*blursize+l]*(g[(-1+height-abs(height-1-(i-shift+k)))*width + width-1-abs(width-1-(j-shift+l))]);
                    blu_tmp+=blur[k*blursize+l]*(b[(-1+height-abs(height-1-(i-shift+k)))*width + width-1-abs(width-1-(j-shift+l))]);
                    
                }
            }
            r2[i*width+j]=(red_tmp);
            g2[i*width+j]=(green_tmp);
            b2[i*width+j]=(blu_tmp);
            red_tmp=0;
            green_tmp=0;
            blu_tmp=0;
    }
    __syncthreads();
    if(i>=0 && i<height && j>=width-shift && j<width)//EDGE RIGHT
    {
            for (int k=0; k<blursize; k++)
            {
                for (int l=0; l<blursize; l++)
                {
                    red_tmp+=blur[k*blursize+l]*(r[(-1+height-abs(height-1-(i-shift+k)))*width + width-1-abs(width-1-(j-shift+l))]);
                    green_tmp+=blur[k*blursize+l]*(g[(-1+height-abs(height-1-(i-shift+k)))*width + width-1-abs(width-1-(j-shift+l))]);
                    blu_tmp+=blur[k*blursize+l]*(b[(-1+height-abs(height-1-(i-shift+k)))*width + width-1-abs(width-1-(j-shift+l))]);
                    
                }
            }
            r2[i*width+j]=(red_tmp);
            g2[i*width+j]=(green_tmp);
            b2[i*width+j]=(blu_tmp);
            red_tmp=0;
            green_tmp=0;
            blu_tmp=0;
    }  */ 
        //PASS TO rgb
    if(i>=0&& i<height && j>=0 && j<=width)
    {
            r[i*width+j]=r2[i*width+j];
            g[i*width+j]=g2[i*width+j];
            b[i*width+j]=b2[i*width+j];
    }
}

/* __global__ void fedge(float* r, float* g, float* b, float* r2, float* g2, float* b2, int width, int height, float*  blur, int blursize) //poor attempt at trying to gain edge filtering back
{
    long long int idx= threadIdx.x + blockIdx.x * blockDim.x;
    int j=int(idx%width);
    int i=int(idx/width);

    float red_tmp=0;
    float green_tmp=0;
    float blu_tmp=0;

    int shift=(blursize-1)/2;
    //EDG TOP
    __syncthreads();
    if(i>=0 && i<shift && j>=0 && j<width)
    {
            for (int k=0; k<blursize; k++)
            {
                for (int l=0; l<blursize; l++)
                {
                    red_tmp+=blur[k*blursize+l]*r[abs(((i-shift)+k))*width + abs((j-shift+l))];
                    green_tmp+=blur[k*blursize+l]*g[abs(((i-shift)+k))*width + abs((j-shift+l))];
                    blu_tmp+=blur[k*blursize+l]*b[abs(((i-shift)+k))*width + abs((j-shift+l))];
                    
                }
            }
            r2[i*width+j]=(red_tmp);
            g2[i*width+j]=(green_tmp);
            b2[i*width+j]=(blu_tmp);
            red_tmp=0;
            green_tmp=0;
            blu_tmp=0;
    }
    __syncthreads();
    //EDG left
    if(i>=0 && i<height && j>=0 && j<shift)
    {
            for (int k=0; k<blursize; k++)
            {
                for (int l=0; l<blursize; l++)
                {
                    red_tmp+=blur[k*blursize+l]*r[abs(((i-shift)+k))*width + abs((j-shift+l))];
                    green_tmp+=blur[k*blursize+l]*g[abs(((i-shift)+k))*width + abs((j-shift+l))];
                    blu_tmp+=blur[k*blursize+l]*b[abs(((i-shift)+k))*width + abs((j-shift+l))];
                    
                }
            }
            r2[i*width+j]=(red_tmp);
            g2[i*width+j]=(green_tmp);
            b2[i*width+j]=(blu_tmp);
            red_tmp=0;
            green_tmp=0;
            blu_tmp=0;
    }
    __syncthreads();
    //EDGE DOWN
    if(i>=height-shift && i<height && j>=0 && j<width)
    {
            for (int k=0; k<blursize; k++)
            {
                for (int l=0; l<blursize; l++)
                {
                    red_tmp+=blur[k*blursize+l]*(r[(-1+height-abs(height-1-(i-shift+k)))*width + width-1-abs(width-1-(j-shift+l))]);
                    green_tmp+=blur[k*blursize+l]*(g[(-1+height-abs(height-1-(i-shift+k)))*width + width-1-abs(width-1-(j-shift+l))]);
                    blu_tmp+=blur[k*blursize+l]*(b[(-1+height-abs(height-1-(i-shift+k)))*width + width-1-abs(width-1-(j-shift+l))]);
                    
                }
            }
            r2[i*width+j]=(red_tmp);
            g2[i*width+j]=(green_tmp);
            b2[i*width+j]=(blu_tmp);
            red_tmp=0;
            green_tmp=0;
            blu_tmp=0;
    }
    __syncthreads();
    if(i>=0 && i<height && j>=width-shift && j<width)//EDGE RIGHT
    {
            for (int k=0; k<blursize; k++)
            {
                for (int l=0; l<blursize; l++)
                {
                    red_tmp+=blur[k*blursize+l]*(r[(-1+height-abs(height-1-(i-shift+k)))*width + width-1-abs(width-1-(j-shift+l))]);
                    green_tmp+=blur[k*blursize+l]*(g[(-1+height-abs(height-1-(i-shift+k)))*width + width-1-abs(width-1-(j-shift+l))]);
                    blu_tmp+=blur[k*blursize+l]*(b[(-1+height-abs(height-1-(i-shift+k)))*width + width-1-abs(width-1-(j-shift+l))]);
                    
                }
            }
            r2[i*width+j]=(red_tmp);
            g2[i*width+j]=(green_tmp);
            b2[i*width+j]=(blu_tmp);
            red_tmp=0;
            green_tmp=0;
            blu_tmp=0;
    } 
    
    __syncthreads();
    //PASS TO rgb
    if(i<shift || i>=height-shift || j<shift || j>=width-shift)
    { 
            r[i*width+j]=r2[i*width+j];
            g[i*width+j]=g2[i*width+j];
            b[i*width+j]=b2[i*width+j];
    }
}*/


__global__ void sepfilterH(float* r, float* g, float* b, float* rnew, float* gnew, float* bnew, int width, int height, float*  blurh, int blursize)
{   __syncthreads();
    long long int idx= threadIdx.x + blockIdx.x * blockDim.x;
    int j=int(idx%width);
    int i=int(idx/width);
    
    float red_tmp=0;
    float green_tmp=0;
    float blu_tmp=0;
    __syncthreads();
    int shift=(blursize-1)/2;
    //MAIN BLUR HORIZONTAL, NO LEFT & RIGHT EDGES
    if(i>=0 && i<height && j>=shift && j<(width-shift))
    {
            for (int l=0; l<blursize; l++)
            {
                red_tmp+=blurh[l]*(rnew[i*width + j-shift+l]);
                green_tmp+=blurh[l]*(gnew[i*width + j-shift+l]);
                blu_tmp+=blurh[l]*(bnew[i*width + j-shift+l]);
          
                //blu_tmp+=blur[k*blursize+l]*float(b[i-(blursize-1)/2+k][j-(blursize-1)/2+l]);
            }
       
            r[i*width+j]=red_tmp;
            g[i*width+j]=green_tmp;
            b[i*width+j]=blu_tmp;

            red_tmp=0;
            green_tmp=0;
            blu_tmp=0;
    }
    __syncthreads();
    //EDGE LEFT HORIZONTAL PART
    //VERTICAL PART WAS DONE IN MAIN BLUR VERTICAL
    if(i>=0 && i<height && j>=0 && j<shift)
    {
            for (int l=0; l<blursize; l++)
            {
                red_tmp+=blurh[l]*(rnew[i*width + abs((j-shift+l))]);
                green_tmp+=blurh[l]*(gnew[i*width + abs((j-shift+l))]);
                blu_tmp+=blurh[l]*(bnew[i*width + abs((j-shift+l))]);
            }
            r[i*width+j]=red_tmp;
            g[i*width+j]=green_tmp;
            b[i*width+j]=blu_tmp;
            red_tmp=0;
            green_tmp=0;
            blu_tmp=0;
    }
    __syncthreads();
    //EDG RIGHT
    if(i>=0 && i<height && j>=(width-shift) && j<width)
    {
            for (int l=0; l<blursize; l++)
            {
                red_tmp+=blurh[l]*(rnew[i*width + width-1-abs(width-1-(j-shift+l))]);
                green_tmp+=blurh[l]*(gnew[i*width + width-1-abs(width-1-(j-shift+l))]);
                blu_tmp+=blurh[l]*(bnew[i*width + width-1-abs(width-1-(j-shift+l))]);
            }
            r[i*width+j]=red_tmp;
            g[i*width+j]=green_tmp;
            b[i*width+j]=blu_tmp;
            red_tmp=0;
            green_tmp=0;
            blu_tmp=0;
    }

}


__global__ void sepfilterV(float* r, float* g, float* b, float* rnew, float* gnew, float* bnew, int width, int height, float* blurv, int blursize) //blur function for spearable filters, much speed wow
{
    __syncthreads();
    long long int idx= threadIdx.x + blockIdx.x * blockDim.x;
    int j=int(idx%width);
    int i=int(idx/width);
    
    //int i=blockIdx.x*blockDim.x+threadIdx.x;
    //int j=blockIdx.y*blockDim.y+threadIdx.y; 

    float red_tmp=0;
    float green_tmp=0;
    float blu_tmp=0;

    int shift=(blursize-1)/2;
    __syncthreads();
    //MAIN BLUR VERTICAL 
    if(i>=shift && i<(height-shift) && j>=0 && j<width)
    {
            for (int k=0; k<blursize; k++)
            {

                //red_tmp+=blurh[k]*float(r[(i-((blursize-1)/2)+k)*width + (j-((blursize-1)/2))]);
                //green_tmp+=blurh[k]*float(g[(i-((blursize-1)/2)+k)*width + (j-((blursize-1)/2))]);
                //blu_tmp+=blurh[k]*float(b[(i-((blursize-1)/2)+k)*width + (j-((blursize-1)/2))]);
                //blu_tmp+=blur[k*blursize+l]*float(b[i-(blursize-1)/2+k][j-(blursize-1)/2+l]);
          
                red_tmp+=blurv[k]*(r[(i-shift)*width+k*width + j]);
                green_tmp+=blurv[k]*(g[(i-shift)*width+k*width + j]);
                blu_tmp+=blurv[k]*(b[(i-shift)*width+k*width + j]);
      
                
            }

            rnew[i*width+j]=red_tmp;
            gnew[i*width+j]=green_tmp;
            bnew[i*width+j]=blu_tmp;

            red_tmp=0;
            green_tmp=0;
            blu_tmp=0;

    }
 
    __syncthreads();
    //EDGE UP USING MIRROR & VERTICAL BLUR
    //SOMEBODY ONCE TOLD ME THAT IF STATEMENTS ARE SLOW
    if(i>=0 && i<shift && j>=0 && j<width)
    {
            for (int k=0; k<blursize; k++)
            {
                red_tmp+=blurv[k]*(r[abs(((i-shift)+k))*width + j]);
                green_tmp+=blurv[k]*(g[abs(((i-shift)+k))*width + j]);
                blu_tmp+=blurv[k]*(b[abs(((i-shift)+k))*width + j]);
            }
            rnew[i*width+j]=red_tmp;
            gnew[i*width+j]=green_tmp;
            bnew[i*width+j]=blu_tmp;
            red_tmp=0;
            green_tmp=0;
            blu_tmp=0;
    }
    __syncthreads();
    //EDGE DOWN
    if(i>=(height-shift) && i<height && j>=0 && j<width)
    {
            for (int k=0; k<blursize; k++)
            {
                red_tmp+=blurv[k]*(r[(-1+height-abs(height-1-(i-shift+k)))*width + j]);
                green_tmp+=blurv[k]*(g[(-1+height-abs(height-1-(i-shift+k)))*width + j]);
                blu_tmp+=blurv[k]*(b[(-1+height-abs(height-1-(i-shift+k)))*width + j]);
            }
            rnew[i*width+j]=red_tmp;
            gnew[i*width+j]=green_tmp;
            bnew[i*width+j]=blu_tmp;
            red_tmp=0;
            green_tmp=0;
            blu_tmp=0;
    }
    

}



void clamp(float* r, float* g, float* b, int width, int height) /////CLAMP TO RGB RANGE
{
    for(unsigned int i=0; i<height; i++)
    {
        for(unsigned int j=0; j<width; j++)
        {
            if (r[i*width+j]>255.0)
            {
                r[i*width+j]=255.0;
            } 
            else if (r[i*width+j]<0)
            {
                r[i*width+j]=0;
            }
            
            if (g[i*width+j]>255.0)
            {
                g[i*width+j]=255.0;
            }
            else if (g[i*width+j]<0)
            {
                g[i*width+j]=0;
            }
            
            
            if (b[i*width+j]>255.0)
            {
                b[i*width+j]=255.0;
            }
            else if (b[i*width+j]<0)
            {
                b[i*width+j]=0;
            }
            
            
        }
    }
}





int main()
{
    int deviceId;
    cudaGetDevice(&deviceId);
    std::cout<<"welcome. input (separated by spaces) your source image path, your chosen destination path and desired operation and parameters for it."<<std::endl<<"0: gaussian blur (parameter: standard deviation (px, float))  ///  1: box blur (parameter: size (px, int, odd))"<<std::endl<<"2: circular box blur (bokeh) (parameter: diameter (px, int, odd))  ///  3: sobel edge detection (parameters: vertical detection, horizontal detection (bool))"<<std::endl<<"4: difference of gaussians edge detection (2 parameters: stdev1, stdev2 (float>0))  ///  5:identity transformation (parameter: size (int)) (useless)"<<std::endl<<"6: unsharp mask simple (parameter: strength (float))  ///  7: unseparated (slow) box blur (parameter: size(px, int, odd))"<<std::endl<<"always type in two numbers as parameters, if the process uses only one parameter set the second one to whatever. 0 or less sets parameter to its default value"<<std::endl;
    std::cout<<"example input: '../samples/default.bmp ../output/test.bmp 0 7 0'"<<std::endl<<std::endl;
    std::string source="default.bmp";
    std::string destination="output/test.bmp";
    int width=0;
    int height=0;
    int sw=0;
    float p1=0;
    float p2=0;
    int p3=0;
    //dim3 threadsperblock(8,8,1);
    dim3 threadsperblock(64,1,1);
    int nchannels=3;
    std::cin>>source>>destination>>sw>>p1>>p2;
    double totaltime=0;
    double elapsedtime=wtime();
    double t0=elapsedtime;
    //check if image exists, is valid format and if yes then query values of width, height and channels
    if (stbi_info(source.c_str(), &width, &height, &nchannels)==0)
    {
        std::cout<<"Invalid source. Aborting mission.";
        return 0;
    }

    if(((long long int)(width)*height)>(2147482137*(threadsperblock.x*threadsperblock.x)))//check if not too big for not using grid stride loops. should be okay for up to many many megapixels to a point where its sligthly sus if you have such a big image
    {
        std::cout<<"listen to me carefully. Tomorrow at 9:37 you have a plane to mexico. I will email you the ticket right away. After you walk out of the airport below a red phone booth is a small container, open it using the secret password: 'hajduszoboszlo'. Therein you will find your new ID card, 3000 pesos and keys to an apartment on the opposite side. From now on your name is Juan Pablo Fernandez Maria FC Barcelona Yanush Sergio Vasilii Shevchenko and you are a russian immigrant from Romania. You work in a barber shop 2 km from the airport. Good luck. Forget about your previous life and maintain a low profile, ditch all of your connections, even with the nvidia customer service."<<std::endl;
        return 0;
    }
    
    //dim3 numberofblocks(((width-1) / threadsperblock.x) + 1, ((height-1) / threadsperblock.y) + 1, 1);
    dim3 numberofblocks((((width*height)-1) / threadsperblock.x) + 1, 1, 1);

    unsigned char *data=stbi_load(source.c_str(), &width, &height, &nchannels, 0);
    //std::cout<<wtime()<<std::endl;    


    float* r;
    float* g;
    float* b;
    float* rnew;
    float* gnew;
    float* bnew; 
    cudaMallocManaged(&r, width*height*sizeof(float));
    cudaMallocManaged(&g, width*height*sizeof(float));
    cudaMallocManaged(&b, width*height*sizeof(float));
    cudaMallocManaged(&rnew, width*height*sizeof(float));
    cudaMallocManaged(&gnew, width*height*sizeof(float));
    cudaMallocManaged(&bnew, width*height*sizeof(float));
    cudaMemPrefetchAsync(r, width*height*sizeof(float), deviceId);
    cudaMemPrefetchAsync(g, width*height*sizeof(float), deviceId);
    cudaMemPrefetchAsync(b, width*height*sizeof(float), deviceId);
    cudaMemPrefetchAsync(rnew, width*height*sizeof(float), deviceId);
    cudaMemPrefetchAsync(gnew, width*height*sizeof(float), deviceId);
    cudaMemPrefetchAsync(bnew, width*height*sizeof(float), deviceId);
   


    //std::cout<<"aeeeeee"<<r[0]<<std::endl;

    for(int i=0; i<height; i++)
    {
        for(int j=0; j<width; j++)
        {
            r[i*width+j]=float(data[i*width*nchannels+j*nchannels]);
            g[i*width+j]=float(data[i*width*nchannels+j*nchannels+1]);
            b[i*width+j]=float(data[i*width*nchannels+j*nchannels+2]);
            //std::cout<<int(data[i*width*nchannels+j*nchannels])<<int(data[i*width*nchannels+j*nchannels+1])<<int(data[i*width*nchannels+j*nchannels+2])<<std::endl;
        }
    }  

    //std::cout<<"aeeeeee"<<std::endl;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////IMYDŻ LOADED
 
    int blursize=3;//keep this odd
    float identity3[]={0,0,0,0,1,0,0,0,0};

switch (sw)
{
case 0:// gausjan
{    ////////GAUS BLUR/////////////    
    float stdev=2;//default
    if (p1>0)
    {
        stdev=p1;
    }
    int gaussize=(2*int(stdev)-1)*3;// 3 sigma wide blur
    if (int(stdev)==0)
    {
        gaussize=3;
    }
    float* gausblur;
    cudaMallocManaged(&gausblur, gaussize*sizeof(float));
    float norm=0;
    for (int i=0;i<gaussize;i++)
    {
        gausblur[i]=1/sqrt(2*3.14159)*exp(-(float)pow(i-(gaussize-1)/2,2)/(2*stdev*stdev));
        norm+=gausblur[i];
    }
    for (int i=0;i<gaussize;i++)
    {
        gausblur[i]/=norm;
        //std::cout<<gausblur[i]<<' ';
    }
    cudaMemPrefetchAsync(gausblur, gaussize*sizeof(float), deviceId);  
    cudaDeviceSynchronize();
    elapsedtime=wtime()-elapsedtime;
    totaltime+=elapsedtime;
    std::cout<<"time elapsed from starting program to beginning of filtering: "<<elapsedtime<<" s."<<std::endl;
    sepfilterV<<<numberofblocks,threadsperblock>>>(r,g,b,rnew,gnew,bnew,width,height,gausblur,gaussize);
    cudaDeviceSynchronize();
    sepfilterH<<<numberofblocks,threadsperblock>>>(r,g,b,rnew,gnew,bnew,width,height,gausblur,gaussize);
    cudaDeviceSynchronize();
    elapsedtime=wtime()-elapsedtime;
    totaltime+=elapsedtime;
    std::cout<<"time elapsed during filtering: "<<elapsedtime<<" s."<<std::endl;
    cudaFree(gausblur);
}
break;
    
  
case 1://box blur (separable)
{
    if (p1>0)
    {
        blursize=(int)(2*ceil(p1/2)-1);
    }
    float* boxblursep;
    cudaMallocManaged(&boxblursep, blursize*sizeof(float));
    for(int i=0;i<blursize;i++)
    {
        boxblursep[i]=1.0/float(blursize);
        //std::cout<<boxblursep[i]<<", ";
    }
    cudaMemPrefetchAsync(boxblursep, blursize*sizeof(float), deviceId);
    cudaDeviceSynchronize();
    elapsedtime=wtime()-elapsedtime;
    totaltime+=elapsedtime;
    std::cout<<"time elapsed from starting program to beginning of filtering: "<<elapsedtime<<" s."<<std::endl; 
    sepfilterV<<<numberofblocks,threadsperblock>>>(r,g,b,rnew,gnew,bnew,width,height,boxblursep,blursize);
    cudaDeviceSynchronize();
    sepfilterH<<<numberofblocks,threadsperblock>>>(r,g,b,rnew,gnew,bnew,width,height,boxblursep,blursize);
    cudaDeviceSynchronize();
    elapsedtime=wtime()-elapsedtime;
    totaltime+=elapsedtime;
    std::cout<<"time elapsed during  filtering: "<<elapsedtime<<" s."<<std::endl;
    cudaFree(boxblursep);
}
break;


  ///////BOKEH (CIRCULAR BOX BLUR) obviously one might say there are no perfect circles in any image composed of pixels and while this is one hundred percent true one may see that they cant actually see the uncircleness beyond a certain point which I shall now take the advantage of
case 2://bokeh very slow not separed
{
    int diameter=5;/////diameter in pixels keep this odd pls
    if (p1>0)
    {
        diameter=(int)(2*ceil(p1/2)-1);
    }
    int radius=(diameter-1)/2;
    float* bokeh;
    cudaMallocManaged(&bokeh, diameter*diameter*sizeof(float));
    float R=0;
    int sb=0;
    for(int i=0; i<diameter; i++)//initialise
    {
        for(int j=0; j<diameter; j++)
        {
            R=sqrt((j-radius)*(j-radius)+(i-radius)*(i-radius));
            if (R<=radius)
            {
                bokeh[i*diameter+j]=1;
                sb+=1;
            }
            else
            {
                bokeh[i*diameter+j]=0;
            }
        }
    }
    for(int i=0; i<diameter; i++)
    {
        for(int j=0; j<diameter; j++)
        {
            bokeh[i*diameter+j]/=float(sb);
            //std::cout<<bokeh[i*diameter+j]<<" ";
        }
    }
    cudaMemPrefetchAsync(bokeh, diameter*diameter*sizeof(float), deviceId);
    cudaDeviceSynchronize();
    elapsedtime=wtime()-elapsedtime;
    std::cout<<"time elapsed from starting program to beginning of filtering: "<<elapsedtime<<" s."<<std::endl;  
    fblur<<<numberofblocks,threadsperblock>>>(r,g,b,rnew,gnew,bnew,width,height,bokeh,diameter);
    cudaDeviceSynchronize();
    elapsedtime=wtime()-elapsedtime;
    std::cout<<"time elapsed during  filtering: "<<elapsedtime<<" s."<<std::endl;
    cudaFree(bokeh);
}
break;

case 3://sobel
{
    bool ver=(bool)int(p1);
    bool hor=(bool)int(p2);
    float* r2;
    float* g2;
    float* b2;
    cudaMallocManaged(&r2,width*height*sizeof(float));
    cudaMallocManaged(&g2,width*height*sizeof(float));
    cudaMallocManaged(&b2,width*height*sizeof(float));
    for(int i=0; i<height; i++)
    {
        for(int j=0; j<width; j++)
        {
            r2[i*width+j]=r[i*width+j];
            g2[i*width+j]=g[i*width+j];
            b2[i*width+j]=b[i*width+j];
        }
    }  
    float gxv[]={1,2,1};
    float gxh[]={1,0,-1};
    float gyv[]={1,0,-1};
    float gyh[]={1,2,1};
    cudaMemPrefetchAsync(&r2,width*height*sizeof(float),deviceId);
    cudaMemPrefetchAsync(&g2,width*height*sizeof(float),deviceId);
    cudaMemPrefetchAsync(&b2,width*height*sizeof(float),deviceId);
    cudaDeviceSynchronize();
    elapsedtime=wtime()-elapsedtime;
    std::cout<<"time elapsed from starting program to beginning of filtering: "<<elapsedtime<<" s."<<std::endl;   
    sepfilterV<<<numberofblocks,threadsperblock>>>(r,g,b,rnew,gnew,bnew,width,height,gxv,3);
    cudaDeviceSynchronize();
    sepfilterH<<<numberofblocks,threadsperblock>>>(r,g,b,rnew,gnew,bnew,width,height,gxh,3);
    cudaDeviceSynchronize();
    sepfilterV<<<numberofblocks,threadsperblock>>>(r2,g2,b2,rnew,gnew,bnew,width,height,gyv,3);
    cudaDeviceSynchronize();
    sepfilterH<<<numberofblocks,threadsperblock>>>(r2,g2,b2,rnew,gnew,bnew,width,height,gyh,3);
    cudaDeviceSynchronize();
    elapsedtime=wtime()-elapsedtime;
    std::cout<<"time elapsed during  filtering: "<<elapsedtime<<" s."<<std::endl;
    
    if(ver && hor)
    {
        for(int i=0; i<height; i++)//combine vertical and horizontal
        {
            for(int j=0; j<width; j++)
            {
                r[i*width+j]=sqrt(r[i*width+j]*r[i*width+j]+r2[i*width+j]*r2[i*width+j]);
                g[i*width+j]=sqrt(g[i*width+j]*g[i*width+j]+g2[i*width+j]*g2[i*width+j]);
                b[i*width+j]=sqrt(b[i*width+j]*b[i*width+j]+b2[i*width+j]*b2[i*width+j]);
            }
        }  
    }
    else if (ver)
    {
        for(int i=0; i<height; i++)//combine vertical and horizontal
        {
            for(int j=0; j<width; j++)
            {
                r[i*width+j]=r2[i*width+j];
                g[i*width+j]=g2[i*width+j];
                b[i*width+j]=b2[i*width+j];
            }
        }
    }
    else if (hor)
    {
        
    }
    else
    {
        for(int i=0; i<height; i++)//combine vertical and horizontal
        {
            for(int j=0; j<width; j++)
            {
                r[i*width+j]=0;
                g[i*width+j]=0;
                b[i*width+j]=0;
            }
        }
    }
    cudaFree(r2);
    cudaFree(g2);
    cudaFree(b2);
}
break;


case 4: //TWO GAUS
{
    float* r2;
    float* g2;
    float* b2;
    float* gausblur;
    float* gausblur2;
    cudaMallocManaged(&r2,width*height*sizeof(float));
    cudaMallocManaged(&g2,width*height*sizeof(float));
    cudaMallocManaged(&b2,width*height*sizeof(float));

    for(int i=0; i<height; i++)
    {
        for(int j=0; j<width; j++)
        {
            r2[i*width+j]=r[i*width+j];
            g2[i*width+j]=g[i*width+j];
            b2[i*width+j]=b[i*width+j];
        }
    }  
    cudaMemPrefetchAsync(&r2,width*height*sizeof(float),deviceId);
    cudaMemPrefetchAsync(&g2,width*height*sizeof(float),deviceId);
    cudaMemPrefetchAsync(&b2,width*height*sizeof(float),deviceId);
    float stdev=2;//default
    if (p1>0)
    {
        stdev=p1;
    }
    int gaussize=(2*int(stdev)-1)*3;// 3 sigma wide blur
    if (int(stdev)==0)
    {
        gaussize=3;
    }
    cudaMallocManaged(&gausblur,gaussize*sizeof(float));
    float norm=0;
    for (int i=0;i<gaussize;i++)
    {
        gausblur[i]=1/sqrt(2*3.14159)*exp(-(float)pow(i-(gaussize-1)/2,2)/(2*stdev*stdev));
        norm+=gausblur[i];
    }
    for (int i=0;i<gaussize;i++)
    {
        gausblur[i]/=norm;
        //std::cout<<gausblur[i]<<' ';
    }
    cudaMemPrefetchAsync(&gausblur,gaussize*sizeof(float),deviceId);
    //second 
    if (p2>0)
    {
        stdev=p2;
    }
    float stdev2=1.2137;
    int gaussize2=(2*int(stdev2)-1)*3;
    if (int(stdev2)==0)
    {
        gaussize2=3;
    }
    cudaMallocManaged(&gausblur2,gaussize2*sizeof(float));
    norm=0;
    for (int i=0;i<gaussize2;i++)
    {
        gausblur2[i]=1/sqrt(2*3.14159)*exp(-(float)pow(i-(gaussize2-1)/2,2)/(2*stdev2*stdev2));
        norm+=gausblur2[i];
    }
    for (int i=0;i<gaussize2;i++)
    {
        gausblur2[i]/=norm;
        //std::cout<<gausblur[i]<<' ';
    }
    cudaMemPrefetchAsync(&gausblur2,gaussize2*sizeof(float),deviceId);

    //difference of gaussians
    cudaDeviceSynchronize();
    elapsedtime=wtime()-elapsedtime;
    std::cout<<"time elapsed from starting program to beginning of filtering: "<<elapsedtime<<" s."<<std::endl;
    cudaDeviceSynchronize();
    sepfilterV<<<numberofblocks,threadsperblock>>>(r,g,b,rnew,gnew,bnew,width,height,gausblur,gaussize);
    cudaDeviceSynchronize();
    sepfilterH<<<numberofblocks,threadsperblock>>>(r,g,b,rnew,gnew,bnew,width,height,gausblur,gaussize);
    cudaDeviceSynchronize();
    sepfilterV<<<numberofblocks,threadsperblock>>>(r2,g2,b2,rnew,gnew,bnew,width,height,gausblur2,gaussize2);
    cudaDeviceSynchronize();
    sepfilterH<<<numberofblocks,threadsperblock>>>(r2,g2,b2,rnew,gnew,bnew,width,height,gausblur2,gaussize2);
    cudaDeviceSynchronize();
    elapsedtime=wtime()-elapsedtime;
    std::cout<<"time elapsed during  filtering: "<<elapsedtime<<" s."<<std::endl;
    for(int i=0; i<height; i++)
    {
        for(int j=0; j<width; j++)
        {
            r[i*width+j]=r[i*width+j]-r2[i*width+j];
            g[i*width+j]=g[i*width+j]-g2[i*width+j];
            b[i*width+j]=b[i*width+j]-b2[i*width+j];
        }
    }
    cudaFree(gausblur); 
    cudaFree(gausblur2);
    cudaFree(r2);
    cudaFree(g2);
    cudaFree(b2);
}
break;


case 5: //identity separable
{
    if (p1>0)
    {
        blursize=(int)(2*ceil(p1/2)-1);
    }
    float* identitysep;
    cudaMallocManaged(&identitysep,blursize*sizeof(float));
    for (int i=0;i<blursize;i++)
    {
        if (i!=(blursize-1)/2)
        {
            identitysep[i]=0.0;
        }
        else
        {
            identitysep[i]=1.0;
        }
    }
    cudaMemPrefetchAsync(&identitysep,blursize*sizeof(float),deviceId);
    cudaDeviceSynchronize();
    elapsedtime=wtime()-elapsedtime;
    std::cout<<"time elapsed from starting program to beginning of filtering: "<<elapsedtime<<" s."<<std::endl;
    cudaDeviceSynchronize();
    sepfilterV<<<numberofblocks,threadsperblock>>>(r,g,b,rnew,gnew,bnew,width,height,identitysep,blursize);
    cudaDeviceSynchronize();
    sepfilterH<<<numberofblocks,threadsperblock>>>(r,g,b,rnew,gnew,bnew,width,height,identitysep,blursize);
    cudaDeviceSynchronize();
    elapsedtime=wtime()-elapsedtime;
    std::cout<<"time elapsed during  filtering: "<<elapsedtime<<" s."<<std::endl;
    free(identitysep);

}
break;


case 6: //unsharp mask 3x3 very slow
{
    float unsharp[9]={0,0,0,0,1,0,0,0,0};
    float crs[]={0,1,0,1,1,1,0,1,0};
    float sharp=6;
    if (p1>0)
    {
        sharp=p1;
    }
    float nor=0;
    for(int i=0; i<3; i++)
    {
        for(int j=0; j<3; j++)
        {
            unsharp[i*3+j]+=(sharp*(identity3[i*3+j]-(crs[i*3+j]/sharp)));
            //std::cout<<unsharp[i*3+j]<<" ";
            nor+=unsharp[i*3+j];
        }
    }
    for(int i=0; i<3; i++)
    {
        for(int j=0; j<3; j++)
        {
            unsharp[i*3+j]/=nor;
           //std::cout<<unsharp[i*3+j]<<" ";
        }
    }
    cudaDeviceSynchronize();
    elapsedtime=wtime()-elapsedtime;
    std::cout<<"time elapsed from starting program to beginning of filtering: "<<elapsedtime<<" s."<<std::endl;
    fblur<<<numberofblocks,threadsperblock>>>(r,g,b,rnew,gnew,bnew,width,height,unsharp,3);
    cudaDeviceSynchronize();
    elapsedtime=wtime()-elapsedtime;
    std::cout<<"time elapsed during  filtering: "<<elapsedtime<<" s."<<std::endl;
}
break;


case 7:// box blur slo
{
    if (p1>0)
    {
        blursize=(int)(2*ceil(p1/2)-1);
    }
    float* boxblur;
    cudaMallocManaged(&boxblur, blursize*blursize*sizeof(float));
    float blooor=float(blursize*blursize);
    for(int i=0;i<blursize;i++)
    {
        for(int j=0;j<blursize;j++)
        {
            boxblur[i*blursize+j]=1.0/blooor;
        }
    }
    cudaDeviceSynchronize();
    elapsedtime=wtime()-elapsedtime;
    std::cout<<"time elapsed from starting program to beginning of filtering: "<<elapsedtime<<" s."<<std::endl;
    fblur<<<numberofblocks,threadsperblock>>>(r,g,b,rnew,gnew,bnew,width,height,boxblur,blursize);
    //cudaDeviceSynchronize();
    //fedge<<<numberofblocks,threadsperblock>>>(r,g,b,rnew,gnew,bnew,width,height,boxblur,blursize);
    cudaDeviceSynchronize();
    elapsedtime=wtime()-elapsedtime;
    std::cout<<"time elapsed during  filtering: "<<elapsedtime<<" s."<<std::endl;    
    cudaFree(boxblur);
}
break;
    
default:
{
    std::cout<<"Unknown option. Aboring program.";
    return 0;
}
break;

}


    

/////////IDENTITY TRANSFORMATION SEPARABLE/////////////////////////////////////// 

    //float* identity3=(float*)malloc(9*sizeof(float));
    //memcpy((int**)id3,identity3,9*sizeof(float));

///////BOX BLUR///////////////////////////////////////////////////////////    


        //std::cout<<boxblursep[1]; 
    /*    std::cout<<std::endl;        
    for(int i=0;i<blursize;i++)
    {
        for(int j=0;j<blursize;j++)
        {
            std::cout<<boxblur[i*blursize+j];
        }

    }*/
    //float boxblur[3][3]={{0,0,0},{0,1,0},{0,0,0}};

 //    std::cout<<"aeeeeee"<<std::endl;




    










//CHOOSE WISELY

//fblur(r,g,b,width,height,boxblur,3);




//fblur(r,g,b,width,height,unsharp,3);









/* 



 //whatever the hell this is
    //std::cout<<"aeeeeee"<<std::endl;

    
    /*for(int i=0; i<height; i++)
    {
        for(int j=0; j<width; j++)
        {
            r[i*width+j]=r[i*width+j]+r2[i*width+j];
            g[i*width+j]=g[i*width+j]+g2[i*width+j];
            b[i*width+j]=b[i*width+j]+b2[i*width+j];
            //std::cout<<int(data[i*width*nchannels+j*nchannels])<<int(data[i*width*nchannels+j*nchannels+1])<<int(data[i*width*nchannels+j*nchannels+2])<<std::endl;
        }
    }  */





    cudaMemPrefetchAsync(r, width*height*sizeof(float), cudaCpuDeviceId);
    cudaMemPrefetchAsync(g, width*height*sizeof(float), cudaCpuDeviceId);
    cudaMemPrefetchAsync(b, width*height*sizeof(float), cudaCpuDeviceId);


cudaDeviceSynchronize();

//////////////////////////////CLAMP/////////
    clamp(r,g,b,width,height);

//////////OUTPUT/////////////////OUTPUT////////////OUTPUT////////////////////////////////////////OUTPUT////////
    for(unsigned int i=0; i<height; i++)
    {
        for(unsigned int j=0; j<width; j++)
        {
            data[i*width*nchannels+j*nchannels]=(int)round(r[i*width+j]);
            data[i*width*nchannels+j*nchannels+1]=(int)round(g[i*width+j]);
            data[i*width*nchannels+j*nchannels+2]=(int)round(b[i*width+j]);
        }
    } 
    
   /*     for(unsigned int i=0; i<height; i++)
    {
        for(unsigned int j=0; j<width; j++)
        {
            data[i*j*nchannels]=r[i][j];
            data[i*j*nchannels+1]=g[i][j];
            data[i*j*nchannels+2]=b[i][j];
        }
    }  */
    //std::cout<<int(data[2]);
    //std::cout<<std::endl;





///////////////////////////WRITE
    stbi_write_bmp(destination.c_str(),width,height,nchannels, data);


/////////F////////////////////REEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

cudaFree(rnew);
cudaFree(gnew);
cudaFree(bnew); 

cudaFree(r);
cudaFree(g);
cudaFree(b);



    cudaDeviceSynchronize();
    elapsedtime=wtime()-elapsedtime;
    totaltime+=elapsedtime;
    //t0=wtime()-t0;
    //std::cout<<"time elapsed during the entire program execution meethod old: "<<t0<<" s."<<std::endl;
    std::cout<<"time elapsed during the last part of program: "<<elapsedtime<<" s."<<std::endl;  
    std::cout<<"time elapsed during the entire program execution: "<<totaltime<<" s."<<std::endl;
    return 0;
}