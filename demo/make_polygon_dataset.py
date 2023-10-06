import numpy as np
import os
import hashlib
import sys
def write_polygon(output_file,rng,w,h):
  nb_vertices = rng.integers(3,11)
  output_file.write(f'{nb_vertices} ')
  alphas = np.sort(rng.uniform(0,np.pi*2,size=nb_vertices))
  d = rng.uniform(min(h,w)/2,max(h,w)/2,size=nb_vertices)
  x = np.cos(alphas)*d
  y = np.sin(alphas)*d
  bx = np.mean(x)
  by = np.mean(y)
  #alphas = np.arctan((y-by)/(x-bx))
  r = ((x-bx)**2+(y-by)**2)**0.5
  assert not np.any(np.isnan(x-bx))
  assert not np.any(np.isnan(y-by))
  assert not np.any(np.isnan(r))
  assert not np.any(np.isnan((x-bx)/r))
  assert not np.any(np.isnan((y-by)/r))
  
  alphas = np.where(y<by,-np.arccos((x-bx)/r),np.arccos((x-bx)/r))
  assert not np.any(np.isnan(alphas)),f'{np.min((x-bx)/r)} {np.max((x-bx)/r)}'
  order = np.argsort(alphas)
  x = x[order]
  y = y[order]
  yshift = y[np.arange(1,y.shape[0]+1)-y.shape[0]]
  xshift = x[np.arange(1,x.shape[0]+1)-x.shape[0]]
  A = 0.5*np.sum((y+yshift)*(x*yshift-xshift*y))
  cx = np.sum((x+xshift)*(x*yshift-xshift*y))/(6*A)
  cy = np.sum((y+yshift)*(y*xshift-yshift*x))/(6*A) 
  x1 = x
  y1 = y
  x= x-cx
  y= y-cy
  #print(alphas,d,bx,by,x1,y1,((x1-bx)**2+(y1-by)**2)**0.5,min(h,w),max(h,w),nb_vertices)
  vertices = []
  for i in range(nb_vertices):
    vertices.append(f'{x[i]} {y[i]}')
  output_file.write(' '.join(vertices))

def write_star_polygon(output_file,rng,w,h):
  s=max(h,w)/2
  nb_branchs = rng.integers(3,7)
  nb_vertices = nb_branchs*2
  output_file.write(f'{nb_vertices} ')
  alphas = np.linspace(0,np.pi*2,1+nb_vertices)[:-1]
  d = np.where(np.arange(nb_vertices)%2==0,np.full((nb_vertices,),s), np.full((nb_vertices,),s/3))
  x = np.cos(alphas)*d
  y = np.sin(alphas)*d
  vertices = []
  for i in range(nb_vertices):
    vertices.append(f'{x[i]} {y[i]}')
  output_file.write(' '.join(vertices))
def write_diamond(output_file,rng,w,h):
  output_file.write(f'0 1 {max(w,h)/2}')
def write_diamond_polygon(output_file,rng,w,h):
  output_file.write(f'4 {-w/2} 0 0 {h/2} {w/2} 0 0 {-h/2}')
def write_rectangle_polygon(output_file,rng,w,h):
  output_file.write(f'4 {-w/2} {-h/2} {-w/2} {h/2} {w/2} {h/2} {w/2} {-h/2}')

def write_clover_polygon(output_file,rng,w,h):
  s = max(h,w)/2
  p1 = []
  res=5
  alpha = np.pi/2
  delta= np.pi*(1/2+30/180)
  a = np.array([np.cos(alpha),np.sin(alpha)])*s/2
  for alpha2 in np.linspace(alpha-delta,alpha+delta,res):
    p1.append(np.array([np.cos(alpha2),np.sin(alpha2)])*s/2+a)
  p1.pop()
  alpha += 2*np.pi/3
  a = np.array([np.cos(alpha),np.sin(alpha)])*s/2
  for alpha2 in np.linspace(alpha-delta,alpha+delta-np.pi/6,res):
    p1.append(np.array([np.cos(alpha2),np.sin(alpha2)])*s/2+a)
  alpha += 2*np.pi/3
  p1.append(np.array([-1/4,-0.9])*s)
  p1.append(np.array([+1/4,-0.9])*s)
  a = np.array([np.cos(alpha),np.sin(alpha)])*s/2
  for alpha2 in np.linspace(alpha-delta+np.pi/6,alpha+delta,res):
    p1.append(np.array([np.cos(alpha2),np.sin(alpha2)])*s/2+a)
  p1.pop()
  output_file.write(f'{len(p1)} ')
  p1 = [f'{p[0]} {p[1]}' for p in p1]
  output_file.write(' '.join(p1))
  
def write_heart_polygon(output_file,rng,w,h):
  s = max(h,w)/2
  p1 = []
  res=5
  alpha = np.pi/4
  delta= np.pi*(1/2)
  a = np.array([np.cos(alpha),np.sin(alpha)])*s/2
  for alpha2 in np.linspace(alpha-delta,alpha+delta,res):
    p1.append(np.array([np.cos(alpha2),np.sin(alpha2)])*s/2+a)
  p1.pop()
  alpha += np.pi/2
  a = np.array([np.cos(alpha),np.sin(alpha)])*s/2
  for alpha2 in np.linspace(alpha-delta,alpha+delta,res):
    p1.append(np.array([np.cos(alpha2),np.sin(alpha2)])*s/2+a)
  p1.append(np.array([0,-1])*s)
  output_file.write(f'{len(p1)} ')
  p1 = [f'{p[0]} {p[1]}' for p in p1]
  output_file.write(' '.join(p1))
  
def write_circle_polygon(output_file,rng,w,h):
  output_file.write(f'0 2 {max(h,w)/2}')

def write_spades_polygon(output_file,rng,w,h):
  s = max(h,w)/2
  p1 = []
  res=5
  alpha = -3*np.pi/4
  delta= np.pi*(1/2)
  a = np.array([np.cos(alpha),np.sin(alpha)])*s*0.6
  for alpha2 in np.linspace(alpha-delta,alpha+delta,res):
    p1.append(np.array([np.cos(alpha2),np.sin(alpha2)])*0.4*s+a)

  alpha += np.pi/2
  p1.append(np.array([-1/4,-0.9])*s)
  p1.append(np.array([+1/4,-0.9])*s)
  a = np.array([np.cos(alpha),np.sin(alpha)])*s*0.6
  for alpha2 in np.linspace(alpha-delta,alpha+delta,res):
    p1.append(np.array([np.cos(alpha2),np.sin(alpha2)])*0.4*s+a)
    
  p1.append(np.array([0,1])*s)

  output_file.write(f'{len(p1)} ')
  p1 = [f'{p[0]} {p[1]}' for p in p1]
  output_file.write(' '.join(p1))
  
"""make one dataset per shape"""
def make_polygon_dataset(input_dir,output_dir):
  filenames = os.listdir(input_dir)
  #funcs = [write_diamond_polygon,write_clover_polygon,write_star_polygon,write_rectangle_polygon,write_heart_polygon,write_spades_polygon,write_circle_polygon]
  funcs = {'diamond':write_diamond_polygon,'clover':write_clover_polygon,'star':write_star_polygon,'rectangle':write_rectangle_polygon,'heart':write_heart_polygon,'spades':write_spades_polygon}
  for shape in funcs.keys():
    os.mkdir(os.path.join(output_dir,shape))
    for filename in filenames:
      m = hashlib.sha256()
      m.update(bytes(os.path.basename(filename),"utf-8"))
      b = m.digest()
      rng = np.random.default_rng([ int.from_bytes(b[4*i:4*(i+1)],"big") for i in range(8)])
      with open(os.path.join(input_dir,filename),'rt') as input_file:
        with open(os.path.join(output_dir,shape,os.path.basename(filename)),'wt') as output_file:
          input_lines = input_file.readlines()
          output_file.write(input_lines[0])
          for line in input_lines[1:]:
            values = line[:-1].split(' ')
            output_file.write(' '.join(values[:2]))
            output_file.write(' ')
            funcs[shape](output_file,rng,float(values[2]),float(values[3]))
            output_file.write(os.linesep)

""" make a datasets where each node can have a different shape""" 
def make_polygon_dataset_multiple_shapes(input_dir,output_dir):
  filenames = os.listdir(input_dir)
  #funcs = [write_diamond_polygon,write_clover_polygon,write_star_polygon,write_rectangle_polygon,write_heart_polygon,write_spades_polygon,write_circle_polygon]
  funcs = {'diamond':write_diamond_polygon,'clover':write_clover_polygon,'star':write_star_polygon,'rectangle':write_rectangle_polygon,'heart':write_heart_polygon,'spades':write_spades_polygon}
  for filename in filenames:
    m = hashlib.sha256()
    m.update(bytes(os.path.basename(filename),"utf-8"))
    b = m.digest()
    rng = np.random.default_rng([ int.from_bytes(b[4*i:4*(i+1)],"big") for i in range(8)])
    with open(os.path.join(input_dir,filename),'rt') as input_file:
      with open(os.path.join(output_dir,shape,os.path.basename(filename)),'wt') as output_file:
        input_lines = input_file.readlines()
        output_file.write(input_lines[0])
        for line in input_lines[1:]:
          values = line[:-1].split(' ')
          output_file.write(' '.join(values[:2]))
          output_file.write(' ')
          rng.choice(funcs,size=1)[0](output_file,rng,float(values[2]),float(values[3]))
          output_file.write(os.linesep)
  
def convert_to_polygon(input_dir,output_dir):
  filenames = os.listdir(input_dir)
  for filename in filenames:
    m = hashlib.sha256()
    m.update(bytes(os.path.basename(filename),"utf-8"))
    b = m.digest()
    rng = np.random.default_rng([ int.from_bytes(b[4*i:4*(i+1)],"big") for i in range(8)])
    with open(os.path.join(input_dir,filename),'rt') as input_file:
      with open(os.path.join(output_dir,os.path.basename(filename)),'wt') as output_file:
        input_lines = input_file.readlines()
        output_file.write(input_lines[0])
        for line in input_lines[1:]:
          values = line[:-1].split(' ')
          output_file.write(' '.join(values[:2]))
          w = float(values[2])
          h = float(values[3])
          output_file.write(' 4 ')
          output_file.write(' '.join(list(map(str,[-w/2,-h/2,w/2,-h/2,w/2,h/2,-w/2,h/2]))))
          output_file.write(os.linesep)

if __name__ == '__main__':
  #for pretier shape approxuimation(heart spades clover) increase variable res
  #make one dataset per shape
  make_polygon_dataset(sys.argv[1],sys.argv[2])
  #use next call to generate one dataset with mulitple shape instead
  #make_polygon_dataset_multiple_shapes(sys.argv[1],sys.argv[2])
