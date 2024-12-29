function beijihanshu=beijihanshu_pow2p(y,rou,u,p,lam)
beijihanshu=(y/rou).^(2-p).*pdf_ed(y,rou,u,p,lam);
end