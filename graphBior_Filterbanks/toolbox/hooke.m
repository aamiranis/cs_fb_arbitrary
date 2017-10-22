function [xopt,Opt,Nav]=hooke(Sx,guess,ip,Lb,Ub,problem,tol,mxit,stp,amp,red,varargin)
%   Unconstrained optimization using Hooke & Jeeves.
%
%   [xopt,Opt,Nav]=hkjeeves(Sx,guess,ip,Lb,Ub,problem,tol,mxit,stp,amp,red,par1, par2,...)
%
%   Sx			: objective function
%   guess		: initial point
%   par			: parameters needed to function
%   ip			: (0): no plot (default), (>0) plot figure ip with pause, (<0) plot figure ip
%   Lb, Ub	: lower and upper bound vectors to plot (default = guess*(1+/-2))
%   problem	: (-1): minimum (default), (1): maximum
%   tol			: tolerance (default = 1e-4)
%   mxit		: maximum number of stages (default = 50*(1+4*~(ip>0)))
%   stp			: stepsize vector for the independent variables (default = max(0.01*abs(guess+~guess),0.1))
%   amp			: stepsize enlargement factor (1,oo) (default = 1.5)
%   red			: stepsize reduction factor (0,1) (default = 0.5)
%   xopt		: optimal point
%   Opt			: optimal value of Sx
%   Nav			: number of objective function evaluations

%   Copyright (c) 2001 by LASIM-DEQUI-UFRGS
%   $Revision: 1.0 $  $Date: 2001/07/05 21:10:15 $
%   Argimiro R. Secchi (arge@enq.ufrgs.br)
%
%
%   Modified by Giovani Tonel (giotonel@enq.ufrgs.br) - April 2007 



%   direction = - 1 --> reverse
%                 1 --> direct
%   idx: variables index vector with enlarged stepsize
%   top: end of idx list
%   bottom: beggin of idx list
%   next: index of enlarged variable


 
 if nargin < 2,
   error('hkjeeves requires two input arguments ''Sx,guess''');
 end
 if nargin < 3 | isempty(ip),
   ip=0;
 end
 if nargin < 4 | isempty(Lb),
   Lb=-guess-~guess;
 end
 if nargin < 5 | isempty(Ub),
   Ub=2*guess+~guess;
 end
 if nargin < 6 | isempty(problem),
   problem=-1;
 end
 if nargin < 7 | isempty(tol),
   tol=1e-4;
 end
 if nargin < 8 | isempty(mxit),
   mxit=50*(1+4*~(ip>0));
 end
 if nargin < 9 | isempty(stp),
   stp=max(0.01*abs(guess+~guess),0.1);
 else
   stp=abs(stp(:));
 end
 if nargin < 10 | isempty(amp) | amp <= 1,
   amp=1.5;
 end
 if nargin < 11 | isempty(red) | red <= 0 | red >= 1,
   red=0.5;
 end

% guess=guess(:);

 yo= feval(Sx,guess, varargin{:} )*problem;
 n=size(guess,1);

  
 x=guess;
 xopt=guess;
 Opt=yo;
 it=0;
 Nav=1;
 top=0;
 bottom=0;
 idx=zeros(n+1,1);
 idx(bottom+1)=top;
 
 while it < mxit,
  next=bottom;
  norma=0;
                     % exploration
  for i=1:n,
   stp_i = stp(i);
   
   for direction=[1 -1],
     x(i)=xopt(i)+stp_i*direction;
     y= feval(Sx,x, varargin{:} )*problem;
     Nav=Nav+1;
     
     if y > yo,     % success
       yo=y;
       if ip & n == 2,
         plot([x(1) xopt(1)],[x(2) xopt(2)],'r');
         if ip > 0,
           disp('Pause: hit any key to continue...'); pause;
         end
       end

       xopt(i)=x(i);
       idx(next+1)=i;
       next=i;
       break;
     end
   end  
  
   x(i)=xopt(i);
   norma=norma+stp_i*stp_i/(x(i)*x(i)+(abs(x(i))<tol));
  end
   
  it=it+1;

  if sqrt(norma) < tol & abs(yo-Opt) < tol*(0.1+abs(Opt)),
    break;
  end
                  % progression
  if next==bottom,
    stp=stp*red;
  else
    good=1;
    idx(next+1)=top;

    while good,
      Opt=yo;
    
      next=idx(bottom+1);
      while next ~= top,
        x(next)=x(next)+amp*(x(next)-guess(next));
        guess(next)=xopt(next);
        next=idx(next+1);
      end
    		  
      if idx(bottom+1) ~= top,
        y= feval(Sx, x, varargin{:} )*problem;
        
        Nav=Nav+1;
      
        if y > yo,
          yo=y;
          good=1;
        
          if ip & n == 2,
            plot([x(1) xopt(1)],[x(2) xopt(2)],'r');
            if ip > 0,
              disp('Pause: hit any key to continue...'); pause;
            end
          end
        else
          good=0;
        end

        next=idx(bottom+1);
        while next ~= top,
          if good,
            xopt(next)=x(next);
          else
            x(next)=xopt(next);
          end
          next=idx(next+1);
        end
      end
    end
  end
 end

 Opt=yo*problem;
 
 if it == mxit,
   disp('Warning Hkjeeves: reached maximum number of stages!');
 elseif ip & n == 2,
   h2=plot(xopt(1),xopt(2),'r*');
   legend([h1,h2],'start point','optimum');
 end
 
disp('Optimization by hkjeeves terminated!')