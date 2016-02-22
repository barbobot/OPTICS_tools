function [OutCellArray, Outconfig]=Map2Zernike(dmap,terms)
%MAP2ZERNIKE computes a set of zernike coefficients from a map.  
% Input parameters are described in ZEMAX order, and then the polynomials
% are recursively generated for fitting.  Both annular and standard
% zernikes may be fit.
%
% Zernikes are only orthogonal over a smooth unit circle which can lead to
% errors when working with discrete data in the form of a matrix. To
% minimize this error, two processing options exist. A weighting map can
% act to smooth the surface [see "Fitting high-order Zernike polynomials to
% finite data", Lewis et al., paper 8493-13, Proceedings of SPIE Volume
% 8493, 2012 ]. Additionally, the process of Gram-Schmidt
% orthonormalization can be applied to force the functions to be
% orthogonal.
%
%Inputs: Map
%Outputs: Zernike
%
%Author: Benjamin Lewis
%Email: blewis@email.arizona.edu
%
% Version history:
% 10/30/11 - Version 2
% 11/14/11 - Version 2.1


%dmap = PME_y;
ztype = 'Standard'; epsilon=0;
        
x=linspace(-1, 1, size(dmap,2));
y=linspace(-1, 1, size(dmap,1));

[x,y]=meshgrid(x,-1*y');

[theta, rho] = cart2pol(x,y);

for ii=1:length(theta(:))
    if theta(ii)>pi
        theta(ii)=theta(ii)-2*pi;
    end
end
weightmap=ones(size(dmap,1),size(dmap,2));
weightmap=WeightMapper(dmap);

theta(isnan(dmap))=[];
rho(isnan(dmap))=[];

dmap(isnan(dmap))=[];

%if isempty(weightmap)
 %   [zvec0,mvec,nvec]=ComputeCoeffsGS(rho, theta, terms, gsterms, epsilon, dmap);
%else
 %   [zvec0,mvec,nvec]=ComputeCoeffs(rho, theta, terms, epsilon, dmap, weightmap);
%end

%OutCellArray{1}.data=[zvec0(:),mvec(:),nvec(:)];
    function wm=WeightMapper(dmap)
         wm=ones(size(dmap,1),size(dmap,2));

         wm(isnan(dmap))=0;

         wm=conv2(wm,[1,1,1;1,0,1;1,1,1],'same');

        if sum(wm(wm==3))<=2
             wm(wm<=2)=0;
             wm(wm<6&wm>2)=.5;
             wm(wm>=6)=1;
         else
             wm(wm<=2)=0;
             wm(wm==3)=.5;
             wm(wm==4)=.5;
             wm(wm==5)=1;
             wm(wm==6)=1;
             wm(wm>=6)=1;
         end

         wm(isnan(dmap))=0;
     end

%    function [wm,rho,theta]=moveweight(ii,jj,wm,dmap,rho,theta)
    %function wm=moveweight(ii,jj,wm,dmap)
        
%         if ii==89&&jj==67
%             pause;
%         end
    
 %       if jj-1>0&&wm(ii,jj-1)~=0&&~isnan(dmap(ii,jj-1)), ...
  %              ln=1; else ln=0; end
        
   %     if jj+1<=size(wm,2)&&wm(ii,jj+1)~=0&&~isnan(dmap(ii,jj+1)), ...
%                rn=1; else rn=0; end
 %       
  %      if ii-1>0&&wm(ii-1,jj)~=0&&~isnan(dmap(ii-1,jj)), ...
   %             un=1; else un=0; end
    %    
%        if ii+1<=size(wm,1)&&wm(ii+1,jj)~=0&&~isnan(dmap(ii+1,jj)), ...
 %               dn=1; else dn=0; end
  %      
   %     nbrs=ln+rn+un+dn;
    %    
     %   if nbrs~=0
      %      w=wm(ii,jj)/nbrs;
       %     
 %           if ln
  %              wm(ii,jj-1)=wm(ii,jj-1)+w;
   %             [a,b]=pol2cart(theta(ii,jj-1),rho(ii,jj-1));
    %            [c,d]=pol2cart(theta(ii,jj),rho(ii,jj));
     %           [tt,tr]=cart2pol((a+w*c)/(1+w),(b+w*d)/(1+w));
      %          rho(ii,jj-1)=tr;
       %         theta(ii,jj-1)=tt;
        %    end
            
        %   if rn
%                wm(ii,jj+1)=wm(ii,jj+1)+w;
 %               [a,b]=pol2cart(theta(ii,jj+1),rho(ii,jj+1));
  %              [c,d]=pol2cart(theta(ii,jj),rho(ii,jj));
   %             [tt,tr]=cart2pol((a+w*c)/(1+w),(b+w*d)/(1+w));
    %            rho(ii,jj+1)=tr;
     %           theta(ii,jj+1)=tt;
      %      end
            
       %     if un
        %        wm(ii-1,jj)=wm(ii-1,jj)+w;
%                [a,b]=pol2cart(theta(ii-1,jj),rho(ii-1,jj));
 %               [c,d]=pol2cart(theta(ii,jj),rho(ii,jj));
  %              [tt,tr]=cart2pol((a+w*c)/(1+w),(b+w*d)/(1+w));
   %             rho(ii-1,jj)=tr;
    %            theta(ii-1,jj)=tt;
     %       end
      %      
   %         if dn
    %            wm(ii+1,jj)=wm(ii+1,jj)+w;
     %           [a,b]=pol2cart(theta(ii+1,jj),rho(ii+1,jj));
      %          [c,d]=pol2cart(theta(ii,jj),rho(ii,jj));
       %         [tt,tr]=cart2pol((a+w*c)/(1+w),(b+w*d)/(1+w));
        %        rho(ii+1,jj)=tr;
         %       theta(ii+1,jj)=tt;
          %  end
      %  end
        
     %   wm(ii,jj)=NaN;
   % end

%     function dmap2=nnavg(ii,jj,dmap,dmap2)
%         if jj-1>0&&~isnan(dmap(ii,jj-1)), ...
%                 ln=1; else ln=0; end
%
%         if jj+1<=size(dmap,2)&&~isnan(dmap(ii,jj+1)), ...
%                 rn=1; else rn=0; end
%
%         if ii-1>0&&~isnan(dmap(ii-1,jj)), ...
%                 un=1; else un=0; end
%
%         if ii+1<=size(dmap,1)&&~isnan(dmap(ii+1,jj)), ...
%                 dn=1; else dn=0; end
%
%         nbrs=ln+rn+un+dn;
%
%         if nbrs~=0
%             v=0;
%
%             if ln, v=v+dmap(ii,jj-1); end
%
%             if rn, v=v+dmap(ii,jj+1); end
%
%             if un, v=v+dmap(ii-1,jj); end
%
%             if dn, v=v+dmap(ii+1,jj); end
%
%             dmap2(ii,jj)=v/nbrs;
%         end
%     end

    function [wm,rho,theta]=WeightMapper2(dmap,rho,theta,pixelspacing,epsilon)
        ps=(sqrt(2)/2)*pixelspacing;
        
        wm=zeros(size(dmap));
        
        for aa=1:length(wm(:))
            if rho(aa)>epsilon-ps&&rho(aa)<epsilon+ps
                [xi,yi]=pol2cart(theta(aa),rho(aa));
                
                x=linspace(-pixelspacing/2,pixelspacing/2,10);
                [x,y]=meshgrid(x,-1*x');
                x=x+xi;
                y=y+yi;
                
                [thetatest,rhotest]=cart2pol(x,y);
                
                rhobool=rhotest;
                
                rhobool(rhotest<=epsilon)=0;
                rhobool(rhotest>epsilon)=1;
                
                
                
                if theta(aa)+pi>.0001 %Removes an ambiguity from near the negative of the polar axis
                    theta(aa)=mean(thetatest(rhobool==1));
                end
                
                rho(aa)=mean(rhotest(rhobool==1));
                
                
                
                wm(aa)=sum(rhobool(:))/length(rhobool(:));
            elseif rho(aa)>epsilon+ps&&rho(aa)<1-ps
                wm(aa)=1;
                
            elseif rho(aa)>1-ps&&rho(aa)<1+ps
                [xi,yi]=pol2cart(theta(aa),rho(aa));
                
                x=linspace(-pixelspacing/2,pixelspacing/2,10);
                [x,y]=meshgrid(x,-1*x');
                x=x+xi;
                y=y+yi;
                
                [thetatest,rhotest]=cart2pol(x,y);
                
                rhobool=rhotest;
                
                rhobool(rhotest<=1)=1;
                rhobool(rhotest>1)=0;
                
                
                
                if theta(aa)+pi>.0001 %Removes an ambiguity from near the negative of the polar axis
                    theta(aa)=mean(thetatest(rhobool==1));
                end
                
                rho(aa)=mean(rhotest(rhobool==1));
                
                
                
                wm(aa)=sum(rhobool(:))/length(rhobool(:));
            end
        end
    end


    function [epsilon] = ObscurationFinder(dmap)
        %UNTITLED Summary of this function goes here
        %   Detailed explanation goes here
        
        bwmap=dmap;
        
        bwmap(isnan(dmap))=0;
        bwmap(~isnan(dmap))=1;
        
        closemap=imfill(bwmap);
        
        holemap=xor(closemap,bwmap);
        
        x=linspace(-1, 1, size(dmap,2));
        y=linspace(-1, 1, size(dmap,1));
        
        [x,y]=meshgrid(x,-1*y');
        
        [~, rho] = cart2pol(x,y);
        
        rho(holemap)=NaN;
        
        rho=rho(:);
        rho(isnan(rho))=[];
        
        rho=unique(rho);
        
        %epsilon=rho(1)*(1-1/size(dmap,1)); %correction for overestimation by a pixel
        epsilon=rho(1);
    end


    function [zout,mout,nout]=ComputeCoeffs(rho, theta, maxTerm, epsilon, dmap, weightmap)
        m=0;
        n=0;
        
        tt=1;
        indices(1,:)=[0,0,0];
        
        while tt<maxTerm
            if m==n||m==-1*n
                n=n+1;
                if mod(n,2)==0, m=0; else m=1; end
            else
                m=m+2;
            end
            if m==0
                indices(tt+1,:)=[m, n, (n-abs(m))/2]; %We can't index to 0, so zernike z0 maps to 1
            else
                if mod(tt,2)~=0 %Strangely enough, ZEMAX uses this to decide the order of terms
                    indices(tt+1,:)=[m, n, (n-abs(m))/2];
                    tt=tt+1;
                    indices(tt+1,:)=[-1*m, n, (n-abs(m))/2];
                else
                    indices(tt+1,:)=[-1*m, n, (n-abs(m))/2];
                    tt=tt+1;
                    indices(tt+1,:)=[m, n, (n-abs(m))/2];
                end
            end
            
            tt=tt+1;
        end
        
        indices=indices(1:maxTerm,:);
        
        if size(indices,1)>1
            Maxs=max(abs(indices));
            MaxM=Maxs(1);
            MaxJ=Maxs(3);
        else
            MaxM=abs(indices(1));
            MaxJ=abs(indices(3));
        end
        
        ThisRow=cell(1,MaxM+MaxJ-1); %This is going the Q Polynomials
        LastRow=cell(1,MaxM+MaxJ-1); %This is going the Q Polynomials
        
        thisH=zeros(1,MaxM+MaxJ-1);
        lastH=zeros(1,MaxM+MaxJ-1);
        
        zout=zeros(maxTerm, 1);
        mout=zeros(maxTerm, 1);
        nout=zeros(maxTerm, 1);
        
        
        %         rho=rho(:);
        %         theta=theta(:);
        %         dmap=dmap(:);
        tdmap=dmap;
        
        rho(end+1)=0;
        theta(end+1)=0;
        
        
        X=(2.*((rho.^2-epsilon^2)./(1-epsilon^2)))-1;
        X(isnan(rho))=NaN;
        
        for mm=0:MaxM
            Qsum=zeros(size(rho));
            for jj=0:MaxM+MaxJ-1-(mm-1)
                if mm==0
                    if jj==0
                        TempQ=ones(size(rho));
                        TempQ(isnan(rho))=NaN;
                        ThisRow{jj+1}=TempQ;
                        thisH(jj+1)=(1-epsilon^2)/(2*(2*jj+1));
                    elseif jj==1
                        ThisRow{jj+1}=X;
                        thisH(jj+1)=(1-epsilon^2)/(2*(2*jj+1));
                    else
                        ThisRow{jj+1}=(2*(jj-1)+1).*X.*ThisRow{jj}-(jj-1).*ThisRow{jj-1};
                        ThisRow{jj+1}=ThisRow{jj+1}./((jj-1)+1);
                        thisH(jj+1)=(1-epsilon^2)/(2*(2*jj+1));
                    end
                else
                    a=2*(2*jj+2*mm-1)*lastH(jj+1);
                    b=(jj+mm)*(1-epsilon^2)*LastRow{jj+1}(end);
                    c=a/b;
                    
                    Qsum=Qsum+LastRow{jj+1}.*LastRow{jj+1}(end)./lastH(jj+1);
                    
                    ThisRow{jj+1}=c.*Qsum;
                    
                    a=2*(2*jj+2*mm-1)*LastRow{jj+2}(end);
                    b=(jj+mm)*(1-epsilon^2)*LastRow{jj+1}(end);
                    
                    thisH(jj+1)=-1*(a/b)*lastH(jj+1);
                end
            end
            
            LastRow=ThisRow;
            lastH=thisH;
            
            rhoM=rho.^mm;
            
            zt=[];
            lvec=[];
            
            for ll=1:size(indices,1)
                if abs(indices(ll,1))==mm
                    m=indices(ll,1);
                    n=indices(ll,2);
                    j=indices(ll,3);
                    
                    scale=(1-epsilon^2)/(2*(2*j+mm+1)*thisH(j+1));
                    scale=scale^.5;
                    scale=scale.*rhoM;
                    
                    if m==0 %Multiply by cos or sin if nessicary, and normalize so RMS=1
                        zernike=sqrt(n+1).*scale.*ThisRow{j+1};
                    elseif m>0
                        zernike=sqrt(2)*sqrt(n+1).*scale.*ThisRow{j+1}.*cos(mm.*theta);
                    else
                        zernike=sqrt(2)*sqrt(n+1).*scale.*ThisRow{j+1}.*sin(mm.*theta);
                    end
                    
                    mout(ll)=m;
                    nout(ll)=n;
                    
                    zernike(end)=[];
                    
                    if isempty(zt)
                        zt{1}=zernike;
                        lvec(1)=ll;
                    else
                        zt{end+1}=zernike;
                        lvec(end+1)=ll;
                    end
                    
                end
            end
            
            [terms,tdmap]=GS(zt,tdmap,weightmap);
            
            for qq=1:length(terms)
                zout(lvec(qq))=terms(qq);
            end
        end
        
        
        function zt=proj(z1,z2,weightmap)
            %zt=(nansum(nansum(z1.*z2.*weightmap))/nansum(nansum(z2.*z2.*weightmap)));
            % remove nansum to avoid unnecessary add-on
            product1 = z1.*z2.*weightmap;
            product2 = z2.*z2.*weightmap;
            numerator = sum(product1(~isnan(product1)));
            denominator = sum(product2(~isnan(product2)));
            zt = numerator/denominator;
        end
        
        %Gram-Schmidt
        function [terms,tdmap]=GS(zc,tdmap,weightmap)
            cmat=zeros(length(zc));
            
            for kk=1:length(zc)
                cmat(kk,kk)=1;
                
                for rr=1:kk-1
                    B=proj(zc{kk},zc{rr},weightmap);
                    
                    zc{kk}=zc{kk}-B*zc{rr};
                    
                    cmat(rr,kk)=-1*B;
                end
                
                tmap=zc{kk}.*tdmap.*weightmap;
                
                a=sum(sum(tmap));
                
                z2=zc{kk}.^2.*weightmap;
                
                b=sum(sum(z2));
                
                z=a/b;
                
                tdmap=tdmap-z.*zc{kk};
                
                terms(kk)=z;
            end
            
            terms=cmat*terms';
        end
        
    end





    function [zout,mout,nout]=ComputeCoeffsGS(rho, theta, maxTerm, gsterms, epsilon, dmap)
        m=0;
        n=0;
        
        tt=1;
        indices(1,:)=[0,0,0];
        
        while tt<maxTerm
            if m==n||m==-1*n
                n=n+1;
                if mod(n,2)==0, m=0; else m=1; end
            else
                m=m+2;
            end
            if m==0
                indices(tt+1,:)=[m, n, (n-abs(m))/2]; %We can't index to 0, so zernike z0 maps to 1
            else
                if mod(tt,2)~=0 %Strangely enough, ZEMAX uses this to decide the order of terms
                    indices(tt+1,:)=[m, n, (n-abs(m))/2];
                    tt=tt+1;
                    indices(tt+1,:)=[-1*m, n, (n-abs(m))/2];
                else
                    indices(tt+1,:)=[-1*m, n, (n-abs(m))/2];
                    tt=tt+1;
                    indices(tt+1,:)=[m, n, (n-abs(m))/2];
                end
            end
            
            tt=tt+1;
        end
        
        indices=indices(1:maxTerm,:);
        
        if size(indices,1)>1
            Maxs=max(abs(indices));
            MaxM=Maxs(1);
            MaxJ=Maxs(3);
        else
            MaxM=abs(indices(1));
            MaxJ=abs(indices(3));
        end
        
        ThisRow=cell(1,MaxM+MaxJ-1); %This is going the Q Polynomials
        LastRow=cell(1,MaxM+MaxJ-1); %This is going the Q Polynomials
        
        thisH=zeros(1,MaxM+MaxJ-1);
        lastH=zeros(1,MaxM+MaxJ-1);
        
        zout=zeros(maxTerm, 1);
        mout=zeros(maxTerm, 1);
        nout=zeros(maxTerm, 1);
        
        
        zc=cell(1,gsterms);
        
        %         rho=rho(:);
        %         theta=theta(:);
        %         dmap=dmap(:);
        tdmap=dmap;
        
        rho(end+1)=0;
        theta(end+1)=0;
        
        
        X=(2.*((rho.^2-epsilon^2)./(1-epsilon^2)))-1;
        X(isnan(rho))=NaN;
        
        for mm=0:MaxM
            Qsum=zeros(size(rho));
            for jj=0:MaxM+MaxJ-1-(mm-1)
                if mm==0
                    if jj==0
                        TempQ=ones(size(rho));
                        TempQ(isnan(rho))=NaN;
                        ThisRow{jj+1}=TempQ;
                        thisH(jj+1)=(1-epsilon^2)/(2*(2*jj+1));
                    elseif jj==1
                        ThisRow{jj+1}=X;
                        thisH(jj+1)=(1-epsilon^2)/(2*(2*jj+1));
                    else
                        ThisRow{jj+1}=(2*(jj-1)+1).*X.*ThisRow{jj}-(jj-1).*ThisRow{jj-1};
                        ThisRow{jj+1}=ThisRow{jj+1}./((jj-1)+1);
                        thisH(jj+1)=(1-epsilon^2)/(2*(2*jj+1));
                    end
                else
                    a=2*(2*jj+2*mm-1)*lastH(jj+1);
                    b=(jj+mm)*(1-epsilon^2)*LastRow{jj+1}(end);
                    c=a/b;
                    
                    Qsum=Qsum+LastRow{jj+1}.*LastRow{jj+1}(end)./lastH(jj+1);
                    
                    ThisRow{jj+1}=c.*Qsum;
                    
                    a=2*(2*jj+2*mm-1)*LastRow{jj+2}(end);
                    b=(jj+mm)*(1-epsilon^2)*LastRow{jj+1}(end);
                    
                    thisH(jj+1)=-1*(a/b)*lastH(jj+1);
                end
            end
            
            LastRow=ThisRow;
            lastH=thisH;
            
            rhoM=rho.^mm;
            
            for ll=1:size(indices,1)
                if abs(indices(ll,1))==mm
                    m=indices(ll,1);
                    n=indices(ll,2);
                    j=indices(ll,3);
                    
                    scale=(1-epsilon^2)/(2*(2*j+mm+1)*thisH(j+1));
                    scale=scale^.5;
                    scale=scale.*rhoM;
                    
                    if m==0 %Multiply by cos or sin if nessicary, and normalize so RMS=1
                        zernike=sqrt(n+1).*scale.*ThisRow{j+1};
                    elseif m>0
                        zernike=sqrt(2)*sqrt(n+1).*scale.*ThisRow{j+1}.*cos(mm.*theta);
                    else
                        zernike=sqrt(2)*sqrt(n+1).*scale.*ThisRow{j+1}.*sin(mm.*theta);
                    end
                    
                    mout(ll)=m;
                    nout(ll)=n;
                    
                    zernike(end)=[];
                    
                    
                    
                    tmap=zernike.*tdmap;
                    %tmap(isnan(tmap))=0;
                    
                    %a=trapz(trapz(tmap));
                    a=sum(sum(tmap));
                    
                    
                    
                    z2=zernike.^2;
                    %z2(isnan(z2))=0;
                    
                    %b=trapz(trapz(z2));
                    b=sum(sum(z2));
                    
                    z=a/b;
                    
                    
                    zout(ll)=z;
                    
                    tdmap=tdmap-z.*zernike;
                    
                    
                    
                    
                    
                    if ll<=gsterms
                        zc{ll}=zernike;
                    end
                end
            end
        end
        
        
        tdmap=dmap;
        
        
        function zt=proj(z1,z2)
            %zt=(nansum(nansum(z1.*z2))/nansum(nansum(z2.*z2)));
            % remove nansum to avoid unnecessary add-on
            product1 = z1.*z2;
            product2 = z2.*z2;
            numerator = sum(product1(~isnan(product1)));
            denominator = sum(product2(~isnan(product2)));
            zt = numerator/denominator;
        end
        
        %Gram-Schmidt
        cmat=zeros(gsterms);
        
        for kk=1:gsterms
            cmat(kk,kk)=1;
            
            for jj=1:kk-1
                B=proj(zc{kk},zc{jj});
                
                zc{kk}=zc{kk}-B*zc{jj};
                
                cmat(jj,kk)=-1*B;
            end
            
            tmap=zc{kk}.*tdmap;
            
            %tmap=zc{kk}.*tdmap;
            
            
            %tmap(isnan(tmap))=0;
            
            %a=trapz(trapz(tmap));
            a=sum(sum(tmap));
            
            
            
            z2=zc{kk}.^2;
            
            %z2=zc{kk}.^2;
            
            
            %z2(isnan(z2))=0;
            
            %b=trapz(trapz(z2));
            b=sum(sum(z2));
            
            z=a/b;
            
            tdmap=tdmap-z.*zc{kk};
            
            tzout(kk)=z;
        end
        
        if gsterms~=0
            tzout2=cmat*tzout';
            
            zout(1:gsterms)=tzout2;
        end
    end



end






