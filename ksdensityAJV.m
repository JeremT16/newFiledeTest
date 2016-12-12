% Copyright (C) 2013 - Michael Baudin
% Copyright (C) 2010 - DIGITEO - Michael Baudin
% Copyright (C) 1993 - 1995 - Anders Holtsberg
%
% This file must be used under the terms of the CeCILL.
% This source file is licensed as described in the file COPYING, which
% you should have received as part of this distribution.  The terms
% are also available at
% http:%www.cecill.info/licences/Licence_CeCILL_V2-en.txt


function [f,xi,u]=ksdensity(varargin)
    % Kernel smoothing density estimate
    %  
    % Calling Sequence
    %   ksdensity(x)
    %   ksdensity(x,u)
    %   ksdensity(x,u,positive)
    %   ksdensity(x,u,positive,akernel)
    %   ksdensity(x,u,positive,akernel,gridsize)
    %   f=ksdensity(...)
    %   [f,xi]=ksdensity(...)
    %   [f,xi,u]=ksdensity(...)
    %             
    % Parameters
    % x : a n-by-1 matrix of doubles, the data
    % u : a 1-by-1 matrix of doubles, positive, the width of the kernel (default=Silverman). On output, u is the kernel width used in the plot.
    % positive : a 1-by-1 matrix of doubles, integer value, set to 1 if the outcomes are positive, set to 0 if negative values of X are possible (default=0). Available values are positive=0, 1.
    % akernel : a 1-by-1 matrix of doubles, integer value, the kernel to use (default=1). Available kernels are 1 (Gaussian), 2 (Epanechnikov), 3 (Biweight), 4 (Triangular)
    % gridsize : a 1-by-1 matrix of doubles, positive, the number of points in the density estimate (default=100)
    % f : a m-by-1 matrix of doubles, the density estimate
    % xi : a m-by-1 matrix of doubles, the corresponding values of X
    %
    % Description
    % Compute a kernel density estimate of the data. 
    %
    % The available kernels are:
    %
    % <itemizedlist>
    % <listitem><para>
    % akernel=1 : Gaussian (the default)
    % </para></listitem>
    % <listitem><para>
    % akernel=2 : Epanechnikov
    % </para></listitem>
    % <listitem><para>
    % akernel=3 : Biweight
    % </para></listitem>
    % <listitem><para>
    % akernel=4 : Triangular
    % </para></listitem>
    % </itemizedlist>
    %
    % The Silverman rule for the width u of the kernel is:
    %
    % <latex>
    % u=1.06\hat{\sigma} n^{-1/5},
    % </latex>
    %
    % where <latex>\hat{\sigma}</latex> is the empirical standard 
    % deviation, and n is the size of the sample.
    %
    % On output, <literal>xi</literal> and <literal>f</literal> contain 
    % the kernel density estimate of the data. 
    % For i=1,2,...,gridsize, f(i) is the estimate 
    % of the probability density at xi(i).
    %
    % Examples
    % X=distfun_normrnd(0,1,1000,1);
    % [f,xi,u]=ksdensity(X);
    % gh=scf();
    % histo(X,[],[],1);
    % plot(xi,f,"r-");
    % xtitle("Kernel density estimate","X","Density");
    % legend(["Data","PDF estimate"]);
    %
    % % Set the kernel width
    % X=distfun_normrnd(0,1,1000,1);
    % [f,xi,u]=ksdensity(X,0.5);
    % scf();
    % plot(xi,f,"r-");
    %
    % % Set the number of points
    % X=distfun_normrnd(0,1,1000,1);
    % [f,xi,u]=ksdensity(X,[],[],[],500);
    % scf();
    % histo(X,[],[],1);
    % plot(xi,f,"r-");
    %
    % % Set the kernel
    % scf();
    % X=distfun_normrnd(0,1,1000,1);
    % %
    % subplot(2,2,1);
    % histo(X,[],[],1);
    % [f,xi,u]=ksdensity(X,[],[],1);
    % plot(xi,f,"r-");
    % xtitle("Gaussian Density Estimate","X","Density")
    % legend(["Data","PDF estimate"]);
    % %
    % subplot(2,2,2);
    % histo(X,[],[],1);
    % [f,xi,u]=ksdensity(X,[],[],2);
    % plot(xi,f,"r-");
    % xtitle("Epanechnikov Density Estimate","X","Density")
    % legend(["Data","PDF estimate"]);
    % %
    % subplot(2,2,3);
    % histo(X,[],[],1);
    % [f,xi,u]=ksdensity(X,[],[],3);
    % plot(xi,f,"r-");
    % xtitle("Biweight Density Estimate","X","Density")
    % legend(["Data","PDF estimate"]);
    % %
    % subplot(2,2,4);
    % histo(X,[],[],1);
    % [f,xi,u]=ksdensity(X,[],[],4);
    % plot(xi,f,"r-");
    % xtitle("Triangular Density Estimate","X","Density")
    % legend(["Data","PDF estimate"]);
    %
    % Authors
    % Copyright (C) 2013 - Michael Baudin
    % Copyright (C) 2010 - DIGITEO - Michael Baudin
    % Copyright (C) 1993 - 1995 - Anders Holtsberg

    %[lhs,rhs]=argn()
    lhs=nargin();
    rhs=nargout ();
%     apifun_checkrhs ( 'ksdensity' , rhs , 1:5 )
%     apifun_checklhs ( 'ksdensity' , lhs , 0:3 )
    f=[];
    xi=[];

    %
    %%%%%je pense qu'ici je verrai le premier commit de la mort et qu'un
    %%%%%jour je serai préseident de la république meme si c'est assez dur
    %%%%%car la vie est simplement plus dure qu'un simple blues sur le pont
    %%%%%des arts avec un max de meufs en train de me regardzer jouer de la
    %%%%%guitare en chantant le blues du soir
    x=varargin(1);
    varlist
    x = x(:);
    n = size(x,'*');
    %
    udefault=1.0600000000000001*st_deviation(x)*(n^(-1/5));
    u=apifun_argindefault(varargin,2,udefault)
    positive=apifun_argindefault(varargin,3,0)
    akernel=apifun_argindefault(varargin,4,1)
    gridsize=apifun_argindefault(varargin,5,100)

    if (x<0) 
        error('There is a negative element in X');
    end

    mn1 = min(x);
    mx1 = max(x);
    mn = mn1-(mx1-mn1)/3;
    mx = mx1+(mx1-mn1)/3;

    xi = linspace(mn,mx,gridsize)';
    d = xi(2)-xi(1);
    xh = zeros(xi);
    xa = (x-mn)/(mx-mn)*gridsize;
    for i = 1:n
        il = floor(xa(i));
        a = xa(i)-il;
        xh(il+[1,2]) = xh(il+[1,2])+[1-a,a]';
    end

    % --- Compute -------------------------------------------------

    xk = ((-gridsize:gridsize-1)')*d;
    if akernel==1 
        K = exp(-0.5*(xk/u).^2);
    
    
    elseif akernel==2 
        K = max(0,1-(xk/u).^2/5);
    
    elseif akernel==3 
        c = sqrt(1/7);
        K = (1-(xk/u)*c.^2).^2 .* (bool2s((1-abs(xk/u*c))>0));
    elseif akernel==4 
        c = sqrt(1/6);
        K = max(0,1-abs(xk/u*c));
    end
    K = K/(sum(K,'m')*d*n);
    %v=size(xh)
    f = fft(fft(fftshift(K),-1) .* fft([xh;zeros(v(1),v(2))],-1),1);
    f = real(f(1:gridsize));

    if positive 
        m = sum(bool2s(xi<0));
        f(m+(1:m)) = f(m+(1:m))+f(m:-1:1);
        f(1:m) = 0
        xi(m+[0,1]) = [0,0];
    end
end 