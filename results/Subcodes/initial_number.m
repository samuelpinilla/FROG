function tot_initial = initial_number(N)

if N <= 64
    
    tot_initial = [4, 12, 20];
    
elseif N == 128
    
    tot_initial = [4, 20, 28];

elseif N == 256
    
    tot_initial = [4, 20, 28];
    
elseif N == 512
    
    tot_initial = [4, 20, 36];
    
elseif N == 1024
    
    tot_initial = [4, 24, 44];
    
elseif N == 2048
   
     tot_initial = [4, 32, 48];
     
elseif N == 4096
    
    tot_initial = [8, 32, 60];
    
else
    
    tot_initial = [8, 36, 64];
    
end
    