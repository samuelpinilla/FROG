function retS = SHG_marginal(Mwin, noise, num_spec) 
  
% inputs: Mwin: Frequency marginal of SHG FROG trace

%         noise: if 0, no filtering

%         num_spec: if 1 retuens only one spectra  

    Mw = quickscale(Mwin);
    N = length(Mw);

    if noise
        
        peri = 5;
        Mw = backsub_noise_1d(Mw, peri);

        Mw(real(Mw)<0) = 1e-6; 

        Mw = smooth(Mw,'sgolay');
  
        Mw = super_gaussian_1d(Mw, .9,20);

        Mw(real(Mw)<0) = 1e-6; 
        
    end
    
    
    coeff_rmv = .2;
  
    St2 = ifftc(1/2/pi*Mw);
 
    St = sqrt(St2);

     
    imagSt = imag(St);          realSt = real(St);      % orginal real and imaginary parts of s(t)
   
    imagStc = imagSt;           realStc = realSt;       % corrected real and imaginary parts of s(t)
    
    pos_neg = ones(size(imagSt));                       % postive or negative coefficent for each compponent of St
      
    
    rmsl = [];   rmsr = [];   rms_ex = [];

    %%% %%% %%%
     
    for i_step =  5: N

         dum_realStc = realStc(1:i_step);
         dum_imagStc = imagStc(1:i_step);
 
         dum_realStc(i_step) = -realStc(i_step);
         dum_imagStc(i_step) = -imagStc(i_step);
         
         
         D3realStc = diff(diff(diff(realStc(1:i_step))));
         D3imagStc = diff(diff(diff(imagStc(1:i_step))));
 
 
         dum_D3realStc = diff(diff(diff(dum_realStc)));
         dum_D3imagStc = diff(diff(diff(dum_imagStc)));
         
         D12r = realStc(i_step) - realStc(i_step-1);
         D01r = realStc(i_step-1) - realStc(i_step-2);
         
         dum_D12r = dum_realStc(i_step) - realStc(i_step-1);
         dum_D01r = realStc(i_step-1) - realStc(i_step-2);
         
         D12im = imagStc(i_step) - imagStc(i_step-1);
         D01im = imagStc(i_step-1) - imagStc(i_step-2);
         
         dum_D12im = dum_imagStc(i_step) - imagStc(i_step-1);
         dum_D01im = imagStc(i_step-1) - imagStc(i_step-2);
         
 
         c = [0.09 .425 1];
         
         wp(1) =  abs(D12r)^2 * c(1) + ...
                  abs(D12r - D01r)^2 * c(2)  + ...
                  abs(D3realStc(i_step-3))^2 * c(3)  ;
         
         wn(1) =  abs(dum_D12r)^2 * c(1) + ...
                  abs(dum_D12r - dum_D01r)^2 * c(2) + ...
                  abs(dum_D3realStc(i_step-3))^2 * c(3) ; 
                   
         wp(2) =  abs(D12im)^2 * c(1) + ...
                  abs(D12im-D01im)^2 * c(2) + ...
                  abs(D3imagStc(i_step-3))^2 * c(3) ;
         
         wn(2) =  abs(dum_D12im).^2 * c(1) + ...
                  abs(dum_D12im - dum_D01im)^2 * c(2) + ...
                  abs(dum_D3imagStc(i_step-3))^2 * c(3) ;
 

         if sum(wp) > sum(wn)
                realStc(i_step:end) = -realStc(i_step:end);
                imagStc(i_step:end) = -imagStc(i_step:end);
                pos_neg(i_step:end) = -pos_neg(i_step:end);        
         end
 
    end
    
    [Swc(1,:), rms0(1), rms0p(1)] = left_side(Mw, St, pos_neg);

    [Swc(2,:), rms0(2), rms0p(2)] = right_side(Mw, St, pos_neg);

% % % EXTRA % % % EXTRA % % % EXTRA % % % EXTRA % % % EXTRA % % % EXTRA % % % EXTRA         
 
    % finding local minima and maxima
    
    real_St_min_finder = realSt;
%     real_St_min_finder(1: round(N * coeff_rmv)-1) = 0;       % to reduce the perfomance time 
%     real_St_min_finder(N-round(N * coeff_rmv)+1 :end) = 0;

    [~, minima_ind] = findpeaks(-real_St_min_finder);
    [~, maxima_ind] = findpeaks(real_St_min_finder); 
    % % 
    
    
    minima_ind(minima_ind < round(N * coeff_rmv)) = [];
    minima_ind(minima_ind > N - round(N * coeff_rmv)) = [];
 
    
    maxima_ind(maxima_ind < round(N * coeff_rmv)) = [];
    maxima_ind(maxima_ind > N - round(N * coeff_rmv)) = [];
    
    nix = length(maxima_ind);
    nim = length(minima_ind);
    
    % making them in the maxima and minima in the same order 
    [~, indp1] =  max(abs(realSt));
    locus = find( minima_ind <  indp1 );
   

    % ordering the minima peakes based on the maxima intensity to their
    % right/ not a necessary step
    if length(locus) > 0
        
        differ = locus(end) - find(maxima_ind == indp1);
    
        minima_ind = circshift(minima_ind, -differ);   
    else
        differ = 0;
    end
    
    if differ < 0
        minima_ind(1:abs(differ)) = [];
        maxima_ind(1:abs(differ)) = [];
    elseif differ > 0
        minima_ind(end-abs(differ):end) = [];
        maxima_ind(end-abs(differ):end) = [];
    end
    % %
    
   if nix > nim
        maxima_ind(length(minima_ind)+1:end)=[];
   else
        minima_ind(length(maxima_ind)+1:end)=[];
   end
 
 
    % finding the points with larg imaginary parts corresponding to the
    % minima points in real part
    

    fakeImag1 = abs(imagStc(minima_ind)); 
    fakeImag2 = abs(imagStc(minima_ind-1));
    fakeImag3 = abs(imagStc(minima_ind+1));
 
    remover = fakeImag1 > .35 * max(abs(imagStc))  & fakeImag2 > .35 * max(abs(imagStc)) & fakeImag3 > .35 * max(abs(imagStc));% &...
      
    % removing the points with large imaginary parts 
    minima_ind( remover ) = [];
    maxima_ind( remover ) = [];
    % %
   
    
    fR = abs(realSt);
    fR(indp1+1:end) = 0;
    ind(1) = indp1; 
    pos_neg1 = pos_neg;     pos_neg1(N/2+2:end) = fliplr(pos_neg1(2:N/2));
    pos_neg2 = pos_neg;     pos_neg2(2:N/2)= fliplr(pos_neg2(N/2+2:end));

    
    for i_peak = 1:sum(maxima_ind<N/2+2)
        
        [~, indmax(i_peak)] = max(abs(fR(maxima_ind)));
    
        min_ind_1 = minima_ind( minima_ind < maxima_ind(indmax(i_peak)));
    
        if length(min_ind_1) > 0
            ind(i_peak) = min_ind_1(end);
        elseif i_peak>1
            ind(i_peak) = ind(i_peak-1);
        end

        fR(maxima_ind(indmax(i_peak))) = 0;
 
        % left side
        pos_neg11(i_peak,:) = pos_neg1;
        pos_neg11(i_peak, ind(i_peak) :end) = - pos_neg1(ind(i_peak):end);   
 
        
        % right side
        indr(i_peak) = N/2+1 +(N/2+1-ind(i_peak)) ;
    
        pos_neg22(i_peak,:) = pos_neg2;        
        pos_neg22(i_peak, indr(i_peak)+1:end) = - pos_neg2(indr(i_peak)+1:end);
     
        %         
        
        [Swcl(i_peak, :),rmsl(i_peak), rmslp(i_peak)] = left_side(Mw, St, pos_neg11(i_peak,:)); 
        %
        [Swcr(i_peak, :),rmsr(i_peak), rmsrp(i_peak)] = right_side(Mw, St, pos_neg22(i_peak,:)); 
 
    end
   
  
    if length(rmsl)>1
               
            [indrms,~] = sort(ind(1:2));
            
            pos_neg_ext = pos_neg1;
            
            A = sort([indrms(1), indrms(2)]);
            pos_neg_ext((A(1))+1: (A(2))-1) = -pos_neg_ext((A(1))+1: (A(2))-1);
     
            [Swc_ex1, rms_ex1, rms_exp1] = left_side(Mw, St, pos_neg_ext);
            
            indrms = sort(indr(1:2));

            A = sort([indrms(1), indrms(2)]);
          
            pos_neg_ext = pos_neg2;
            pos_neg_ext((A(1)): (A(2))) = - pos_neg_ext((A(1)): (A(2)));
            [Swc_ex2, rms_ex2, rms_exp2] = right_side(Mw, St, pos_neg_ext);
           
            if rms_ex1 < rms_ex2
                Swc_ex = Swc_ex1;
                rms_ex = rms_ex1;
                rms_exp = rms_exp1;
            else
                Swc_ex = Swc_ex2;
                rms_ex = rms_ex2;
                rms_exp = rms_exp2;
            end
    end
  
        
	if ~isempty(rmsl)
            if length(rmsl)>1
                rms_total = [rms0, rmsr, rmsl, rms_ex];
                rms_total_p = [rms0p, rmsrp, rmslp, rms_exp];
                S_total = [Swc; Swcr; Swcl; Swc_ex];
            else
                rms_total = [rms0, rmsr, rmsl];
                rms_total_p = [rms0p, rmsrp, rmslp];
                S_total = [Swc; Swcr; Swcl];
            end


            [condf, order] = sort(rms_total);
             num_s = length(order);
             order_o = order;

            
            tot_initial = initial_number(N);

            V = 0;
            if num_s > 2     
                V = 3;
            end    

            num_asg = 1;
            if V > 0
                if num_s < 4/3 * V
                    for ii = 1:V+1
                        Swf(ii,:) = quickscale(S_total(order_o(ii),:));
                    end
                elseif  N < 512 
                    [condf_p, order_p] = sort(rms_total_p);
                  
                    [order_union, io, iop] = intersect(order(1:ceil(3*num_s/4)), order_p(1:ceil(3*num_s/4)));
         
                    S_total_u = S_total(order(io),:);
       
                    [condf, order] = sort(rms_total(order(io)));
                    
                    S_total_p = S_total_u(order,:);

                    Swf(1,:) = S_total_p(1,:);
                    spec_fill(1) = 1;
                    for ii = 1:V                   
                        num_asg = num_asg+1; 
                        delta = zeros(1,length(order));
                        for jj = 1: length(order)
                            if ~sum(jj==spec_fill)

                                for kk = 1:num_asg-1
                                     delta_i = mean(abs(quickscale(S_total_p(jj,:)) - quickscale(Swf(kk,:)))); 
                                     delta(jj) = delta(jj)+ delta_i;
                                end

                            else
                              
                                delta(jj) = 0;
                            end
                        end
                        [~,ind] = sort(delta);
                        Swf(ii+1,:) = quickscale(S_total_p(ind(end),:));        %selecting the most diffirent ones
                        spec_fill(ii+1) = ind(end);     
                    end
                else
                    [condf_p, order_p] = sort(rms_total_p);
                    [order_union, io, iop] = intersect(order(1:round(1*num_s/2)), order_p(1:round(1*num_s/2)));
         
                    S_total_u = S_total(order(io),:);
       
                    [condf, order] = sort(rms_total(order(io)));
                    
                    S_total_p = S_total_u(order,:);

                    Swf(1,:) = S_total_p(1,:);
                    spec_fill(1) = 1;
                    for ii = 1:V                   
                        num_asg = num_asg+1; 
                        delta = zeros(1,length(order));
                        for jj = 1: length(order)
                            if ~sum(jj==spec_fill)

                                for kk = 1:num_asg-1
                                     delta_i = mean(abs(quickscale(S_total_p(jj,:)) - quickscale(Swf(kk,:)))); 
                                     delta(jj) = delta(jj)+ delta_i;
                                end

                            else
                              
                                delta(jj) = 0;
                            end
                        end
                        [~,ind] = sort(delta);
                        Swf(ii+1,:) = quickscale(S_total_p(ind(end),:));        %selecting the most diffirent ones
                        spec_fill(ii+1) = ind(end);     
                    end
                end
                
            end
	else
        condf = min(rms0) ;
        if condf == rms0(1)
            Swf(1,:) = Swc(1,:);
            Swf(2,:) = Swc(2,:);
        else
            Swf(1,:) = Swc(2,:);
            Swf(2,:) = Swc(1,:);
        end
	end   
         

    if num_spec == 1
        retS = quickscale(abs(Swf(1,:)));
    else 
        retS = abs(Swf);
    end
    
end
 
 

