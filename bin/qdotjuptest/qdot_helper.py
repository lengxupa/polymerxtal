# utilities for qdot_simtool
import math

def eval_core_simtime(num_atoms, tb_basis):
    
    num_atoms = num_atoms/1e6;
    coeff1 = {'1':0,'2':0,'4':0,'8':0,'16':0,'32':0}
    coeff2 = {'1':0,'2':0,'4':0,'8':0,'16':0,'32':0}
    coeff3 = {'1':0,'2':0,'4':0,'8':0,'16':0,'32':0}
    coeff4 = {'1':0,'2':0,'4':0,'8':0,'16':0,'32':0}
    time = 0
    cores = 1
    simtime = ""
              
    if (tb_basis == "VBCB" or tb_basis == "em"):

        coeff1['1']=2.2370
        coeff2['1']=-5.7323
        coeff3['1']=108.3358
        coeff4['1']=-0.7924

        coeff1['2']=0.5450
        coeff2['2']=-1.4245
        coeff3['2']=58.8979
        coeff4['2']=-0.2164

        coeff1['4']=0.2991
        coeff2['4']=-0.7930
        coeff3['4']=32.7135
        coeff4['4']=-0.1398

        coeff1['8']=0.1961
        coeff2['8']=-0.5175
        coeff3['8']=19.7550
        coeff4['8']=-0.0874

        coeff1['16']=0.3402
        coeff2['16']=-0.9229
        coeff3['16']=17.8022
        coeff4['16']=-0.1690

    elif (tb_basis == "sp3sstar") :

        coeff1['32']=-6.2108
        coeff2['32']=26.0559
        coeff3['32']=18.1968
        coeff4['32']=0.7644       

    elif (tb_basis == "sp3sstar_SO"):

        coeff1['32']=1.9743 
        coeff2['32']=0.2953
        coeff3['32']=37.9150
        coeff4['32']=-0.1674
        
    elif (tb_basis == "sp3d5sstar"):

        coeff1['32']=1.9743 
        coeff2['32']=0.2953
        coeff3['32']=37.9150
        coeff4['32']=-0.1674
        
    elif (tb_basis == "sp3d5sstar_SO"):
        
        coeff1['32']=7.4564   
        coeff2['32']=-6.7151
        coeff3['32']=126.8985
        coeff4['32']=-2.5647
  
    if (tb_basis == "VBCB" or tb_basis == "em"):
        time1 = coeff1['1']*math.pow(num_atoms,3) + coeff2['1']*math.pow(num_atoms,2) + coeff3['1']*num_atoms + coeff4['1']
        time2 = coeff1['2']*math.pow(num_atoms,3) + coeff2['2']*math.pow(num_atoms,2) + coeff3['2']*num_atoms + coeff4['2']
        time4 = coeff1['4']*math.pow(num_atoms,3) + coeff2['4']*math.pow(num_atoms,2) + coeff3['4']*num_atoms + coeff4['4']
        time8 = coeff1['8']*math.pow(num_atoms,3) + coeff2['8']*math.pow(num_atoms,2) + coeff3['8']*num_atoms + coeff4['8']
        time16 = coeff1['16']*math.pow(num_atoms,3) + coeff2['16']*math.pow(num_atoms,2) + coeff3['16']*num_atoms + coeff4['16']
        if (time1 < 5) :
            time = int(time1)+1
            cores = 1
        elif (time2 < 5) :
            time = int(time2*1.5)
            cores = 2
        elif (time4 < 5):
            time = int(time4*1.5)
            cores = 4
        elif (time8 < 5):
            time = int(time8*1.5)
            cores = 8
        else :
            time = int(time16*1.5)
            cores = 3   
            
    else:      
        time['32'] = coeff1['32']*math.pow(num_atoms,3) + coeff2['32']*math.pow(num_atoms,2) + coeff3['32']*num_atoms + coeff4['32'] 
        time = int(time['32']*1.8 + 10)
        cores = 16
        time = time*2
    
    return time, cores


import os
import sys
import importlib.machinery
loader = importlib.machinery.SourceFileLoader('vtk', '/apps/share64/debian7/anaconda/anaconda3-5.1/lib/python3.6/site-packages/vtk/__init__.py')
vtk = loader.load_module('vtk')

def resample_vtk_serial(filenames, dim=[35,35,35], radius = 0.1, mode='max', exp=None):
    fields = []
    for filename in filenames:
        field = []
        reader = vtk.vtkPolyDataReader()
        path = os.path.join(filename)
        reader.SetFileName(path)
        reader.Update();
        datar = reader.GetOutput()
        n = datar.GetNumberOfPoints()
        origin = [0,0,0]
        end = [0,0,0]
        coord = datar.GetPoint(0)
        origin[0], origin[1], origin[2] = coord[:3]
        end[0], end[1], end[2] = coord[:3]
        for i in range(n):
            coord = datar.GetPoint(i)
            for id, val in enumerate(coord[:3]):
                if origin[id] > val:
                    origin[id] = val
                if end[id] < val:
                    end[id] = val
        points = reader.GetOutput().GetPointData().GetArray(0)

        datar.GetPointData().SetScalars(points)

        splatter = vtk.vtkGaussianSplatter();
        splatter.SetModelBounds(origin[0], end[0], origin[1], end[1], origin[2], end[2])
        if vtk.VTK_MAJOR_VERSION <= 5:
            splatter.SetInput(datar)
        else:
            splatter.SetInputData(datar)
        grid_sample = [dim[0], dim[1], dim[2]]
        splatter.SetSampleDimensions(grid_sample[0],grid_sample[1],grid_sample[2]);
        splatter.SetRadius(radius);
        if mode=='sum':
            splatter.SetAccumulationModeToSum()
        else:
            splatter.SetAccumulationModeToMax()
        if exp is not None:
            splatter.NormalWarpingOff()
            splatter.SetExponentFactor (exp)
        splatter.Update()
        data = splatter.GetOutput()

        points = data.GetPointData().GetScalars()
        n = points.GetNumberOfTuples()
        for i in range(n):
            coord = data.GetPoint(i)
            field.append([coord[0], coord[1], coord[2], points.GetTuple1(i)])
        fields.append(field)
    return fields


def interpolate_vtk_serial(filenames, dim=[30,30,30], radius = 0.2, mode='max', exp=None):
    fields = []
    for filename in filenames:
        field = []
        reader = vtk.vtkPolyDataReader()
        path = os.path.join(filename)

        
        reader.SetFileName(path)
        reader.Update();
        datar = reader.GetOutput()

        center = datar.GetCenter()
        bounds = datar.GetBounds()
        length = datar.GetLength()

        volume = vtk.vtkImageData()
        volume.SetDimensions(dim[0],dim[1],dim[2])
        volume.SetOrigin(bounds[0],bounds[2],bounds[4])
        volume.SetSpacing((bounds[1]-bounds[0])/(dim[0]-1),
                         (bounds[3]-bounds[2])/(dim[1]-1),
                         (bounds[5]-bounds[4])/(dim[2]-1))

        points = datar.GetPoints()
        oneScalars = datar.GetPointData().GetArray(0)

        oneData = vtk.vtkPolyData()
        oneData.SetPoints(points)
        oneData.GetPointData().SetScalars(oneScalars)

        gaussianKernel = vtk.vtkGaussianKernel()
        gaussianKernel.SetRadius(1.)
        gaussianKernel.SetSharpness(2.)

        interpolator = vtk.vtkPointInterpolator()
        interpolator.SetInputData(volume)
        interpolator.SetSourceData(oneData)
        interpolator.SetKernel(gaussianKernel)
        interpolator.SetNullPointsStrategyToClosestPoint()
        interpolator.PassPointArraysOn()
        interpolator.Update()

        data = interpolator.GetOutput()
        delta = volume.GetSpacing();
        for i in range(dim[0]):
            x = bounds[0] + i*delta[0]
            for j in range(dim[1]):
                y = bounds[1] + j*delta[1]
                for k in range(dim[2]):
                    z = bounds[2] + k*delta[2]
                    v = data.GetScalarComponentAsDouble(i,j,k,0)
                    field.append([x,y,z,v])
        fields.append(field)        
    return fields
