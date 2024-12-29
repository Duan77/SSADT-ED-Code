function pdf=pdf_ed(y,t,u,p,lam)
 pdf=sqrt(lam./(2*pi*t^(1-p)*y.^p)).*exp(-lam*t*(((y/t).^(2-p)/(1-p)/(2-p)-y/t*u^(1-p)/(1-p)+u^(2-p)/(2-p))));
end
