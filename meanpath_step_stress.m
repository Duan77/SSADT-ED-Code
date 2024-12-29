clear
theta=[0.124266173663911,1.54458138438586,1.01949946485847,0.484549499028519,1.70641297762755];

u=theta(1);
p=theta(2);
lam=theta(3);
b=theta(4);
a=theta(5);

Y1=[0	3.492614043	5.767917665	6.660752002	7.303833425	7.597952322	8.797196628	8.870447645	9.036681181	10.29547709	10.40679711
0	2.078060159	3.036692689	3.27584302	3.859490063	4.18158148	4.890374388	5.027768974	5.117144256	5.679018224	6.113321018
0	4.094945613	5.588550059	5.974969939	6.003553454	6.292890362	6.453779205	6.587799154	6.911120662	7.331813633	7.46654262
0	3.965079339	4.858639607	5.115135946	6.338866187	6.419415872	6.666775193	7.230684594	7.579932155	8.065293861	8.336738445
0	2.381501693	3.043295825	3.673405207	4.765652506	4.825081626	5.411135549	5.702163308	5.980511039	6.339901203	6.531976249
0	1.951797196	3.093275395	4.320243255	5.09172646	5.328999204	5.694844151	6.068563385	6.392318394	6.563301298	7.059586391
0	3.185180539	3.874173637	4.09195999	4.198588327	4.725688635	4.978048436	5.817698515	6.751236786	6.803861647	6.870059282
0	2.728945029	4.043729559	4.390312289	5.050140348	6.072064885	6.572733098	6.736680309	7.426718118	8.548264981	9.096179292
0	2.36698101	3.046530544	3.188483649	3.348585206	3.52740492	4.241509308	4.716238481	5.221253968	5.673434833	6.445320844
0	2.449712016	4.452265298	5.673684807	7.245311616	8.30463138	8.458908117	8.506737019	8.667964889	8.97964578	9.389680174
0	2.894068655	3.31049552	3.597524925	4.60518217	5.634518847	5.831101209	6.312128047	7.180287495	7.650961322	8.064194878
0	1.725399002	2.564178876	3.844232432	3.906846912	4.480729486	4.822251997	5.108918667	5.414615173	5.428071363	5.699838801
0	2.383199022	3.492526879	4.646581176	5.364948804	5.970373936	7.037669751	8.509482201	9.41289234	9.704775986	10.35888606
0	2.044026089	3.283824595	3.807283321	4.509016739	4.834351643	5.554968448	5.802031417	6.243344764	7.47405528	7.494128702
0	2.269029611	3.197604295	4.284934934	4.69688976	5.572686167	6.632060256	7.696305397	8.302644651	8.387993028	8.855069888
0	2.240684554	2.809803959	3.825466486	4.287669824	4.354914675	4.99886678	5.248085409	5.845666423	7.558937696	7.757111591
0	1.93154055	2.460416113	2.922164277	4.130385889	5.135002407	5.611563391	6.112003927	6.940233571	7.206915216	8.172330718
0	2.121657206	2.32748977	2.859534679	2.964971456	3.197891294	3.894231297	4.76168904	5.453335977	5.898825236	6.767408434
0	3.01339511	4.625082969	5.007382075	5.354810157	5.719707339	5.80728625	5.846592968	6.712814346	6.722291827	7.998142146
0	3.473704013	3.932504439	4.166671064	4.445322817	5.25958907	6.393260089	6.400734525	6.763752396	6.96610024	7.052409372
0	2.723529082	3.702851219	5.015108293	5.861338027	6.284226514	6.615915594	7.232330647	7.309319531	7.674537488	7.766112991
0	2.851111448	3.429256438	4.135560322	4.285648083	4.985385039	5.908990084	6.316351267	7.287016921	7.380109697	7.858220513
0	0.959476107	2.272613749	3.805828994	3.989700245	4.69356574	5.20014545	6.068498769	6.720415999	6.947458226	7.422543917
0	2.650198125	3.743498165	4.227068066	4.686715028	5.403659591	5.874413272	6.602360004	6.685564033	6.988175011	7.921707882
0	3.527875976	4.36843759	4.523781174	4.556645972	5.544022281	6.132100172	6.66156573	6.765053086	6.992415999	7.69105115
0	1.942729949	3.151359284	3.380280757	3.726793209	4.773787251	5.41628696	6.579344137	7.33583952	7.849637866	8.658270659
0	1.695058585	2.39083171	4.292478111	4.706303829	4.772061343	5.031805128	5.130937067	5.152808818	5.514110586	5.84048476
0	1.799958022	2.492425155	3.363109784	3.451941165	4.265562075	4.517101452	4.563832502	4.708854027	5.283351424	6.171253468
0	4.356001968	4.716223806	5.106551625	5.885184279	6.293431601	7.236898089	7.262097568	7.375119234	7.449695038	8.079079087
0	2.387145297	3.207077597	3.298045706	3.340954445	4.962188417	5.540335675	5.587642652	5.931442385	6.489198018	6.641129983
];

Y2=[10.40679711	11.44443376	11.88084575	11.99861244	12.88238946	13.04618792	14.42288206	14.90939493	15.25398588	15.64876966	15.73969367
6.113321018	7.606466587	8.48606913	9.571382839	9.705351008	10.03354018	10.43026579	11.44653631	12.08117237	12.09161642	13.01925531
7.46654262	8.21553913	8.627190458	9.033767613	9.784981038	10.2947436	10.75928891	11.40152651	12.49514754	13.84146866	14.73415038
8.336738445	9.69359654	10.37779075	10.73963533	11.47838557	12.62149927	13.10399866	13.60382456	13.73720429	14.13987369	15.45604155
6.531976249	7.596662612	8.544975071	9.203243167	9.498483126	9.926652264	9.99626399	10.06701244	10.22964798	10.48255808	10.83446939
7.059586391	8.669427667	9.478534678	10.31662853	10.87450703	12.18046324	12.50018443	14.10311113	14.3009622	14.34934542	15.18312939
6.870059282	7.689924987	10.13731427	11.29420865	12.35022199	13.28471166	13.63636577	14.21934227	14.70238876	15.27599198	15.6310119
9.096179292	10.40187649	11.66979211	12.22113036	13.30119605	13.56699158	14.45937503	15.0931173	15.8567141	16.60926729	17.08908935
6.445320844	7.904224922	9.144415405	9.968349781	10.97657021	12.44091939	13.02876873	13.5883104	15.12294927	15.59241202	16.37403281
9.389680174	10.76761265	12.36665979	14.53718935	14.7776715	15.15860272	15.51804772	16.49079862	17.49261465	18.2351874	18.66599055
8.064194878	8.924985879	9.697425499	9.851445976	10.04875675	10.68430489	11.15028018	11.85144958	12.23272525	12.47896572	13.07518047
5.699838801	6.582158364	7.645520207	8.859159852	9.220030917	10.61471606	11.4647157	12.191807	12.66563865	12.97749325	13.80200626
10.35888606	11.42167563	11.7487979	12.00720311	12.94006583	13.60481386	13.6151261	14.04144414	14.28410894	14.5470834	15.19776349
7.494128702	7.900759237	8.94678646	9.803636499	10.60494466	11.63634701	12.17335826	12.33170508	12.80868433	12.98261192	13.50661475
8.855069888	9.248334904	10.09562706	10.21594574	10.6910199	11.29442606	12.13266067	13.47499481	14.36351597	14.57579819	15.08539834
7.757111591	8.993905075	9.638976678	10.03889937	10.57676856	11.69263414	11.91720166	12.4578076	12.72278823	13.0841492	13.68753534
8.172330718	9.319745473	10.02902682	10.81128104	12.44824263	12.76786866	13.90962708	13.91762437	14.15160971	14.74045808	14.87312764
6.767408434	8.113000554	8.211703729	8.590784813	9.353765454	10.95836286	12.83904759	14.1788904	14.94368589	15.8165016	16.91978721
7.998142146	9.502906657	10.53560362	12.13994648	12.49717413	14.17975049	14.60974451	15.02790642	15.47843702	15.54489507	15.70663563
7.052409372	7.848221731	8.063777246	8.948568795	10.09880652	10.38634184	11.02241085	11.71346433	12.36385814	13.544107	13.78845438
7.766112991	8.261718712	9.506925283	12.34437172	14.52702722	15.91891724	16.70344314	16.81787562	17.55493191	17.67628664	18.258099
7.858220513	9.287312831	10.46370327	11.35781924	12.21155758	14.04182871	14.49919774	15.82590542	16.59376251	17.2022948	17.48502182
7.422543917	9.805029281	10.83744906	11.73132245	12.02078722	13.2779663	14.1930623	14.49393867	14.63718439	15.45687449	16.4963768
7.921707882	9.27537324	9.635760126	9.798037475	11.27757545	12.1361648	13.03276507	13.82424327	14.75792574	14.920166	15.75606033
7.69105115	8.717404432	9.899680736	11.00202396	12.43215542	14.13134431	15.35112468	15.77910809	16.18685708	16.79358452	17.68272621
8.658270659	10.28094134	10.9617437	11.32997091	12.25089222	12.46000807	13.34889855	13.9985559	14.37410558	15.02109395	15.17895
5.84048476	5.926722733	6.614744266	8.262409919	10.43895357	11.04550729	11.16680392	12.30670328	12.50970967	13.43364255	13.55090502
6.171253468	6.651955266	7.496745847	8.161033162	8.791502142	9.465980302	10.0268074	10.6421263	11.15292682	11.44623846	12.09987271
8.079079087	10.55275683	11.48507046	12.83698015	13.37679065	14.6529905	15.06873478	15.54570575	16.37218364	17.38189108	18.23450576
6.641129983	7.994136205	9.958449404	11.37196137	12.41544656	13.0568005	13.40554498	14.79786238	15.42789229	15.89982425	16.16710918
];

Y3=[15.73969367	16.54776157	18.22061267	19.73284107	21.19133476	21.6880623	24.19136252	25.28014363	25.29971973	26.25149246	26.95017436
13.01925531	13.57803668	14.54213079	15.98304514	16.71540168	17.36706205	19.00732492	19.42832916	19.8954549	21.09921506	21.94701207
14.73415038	16.24979031	16.45959946	17.60605471	19.14933606	19.89685029	20.41542925	21.49936017	23.05475393	23.35296067	24.64009192
15.45604155	16.27719305	17.32830347	17.53773002	18.73514853	19.6519297	20.63126374	20.65421236	21.07075362	21.15129284	22.18026735
10.83446939	12.0967598	12.98759437	14.77697015	15.40421318	16.36230666	17.05822957	17.98578758	20.59745264	21.7101408	22.62095066
15.18312939	15.94901675	16.80778215	17.67221943	19.93096583	22.65585489	24.34466709	24.46365432	24.59168482	25.33818811	25.73966092
15.6310119	16.29374017	18.07826767	18.91557584	20.85993783	22.19553769	23.17263546	24.77493138	25.51629574	26.56666733	26.86525976
17.08908935	17.2179892	17.93781018	19.61887482	20.13194284	20.45446848	21.26316459	21.67290708	22.48576041	24.64772494	25.09849749
16.37403281	17.4207671	18.19812274	18.60803046	19.98285732	20.7681295	21.59631149	22.27186726	22.63315725	23.00275759	23.26948344
18.66599055	19.43744634	19.81489786	20.50251802	20.87841621	21.26414502	22.75526066	24.60362993	25.90689572	26.35810474	26.55601034
13.07518047	13.75074208	14.23964144	14.31583448	15.23328633	16.02169939	17.11920879	17.32713241	17.81222955	18.93298605	19.56708284
13.80200626	15.11809257	16.83417268	18.08496975	18.45164784	20.05815843	21.91808104	22.34182949	22.74224907	23.90203436	24.62341126
15.19776349	16.0727499	16.08597856	16.89073045	17.45281667	17.86944957	19.76818365	20.76987735	21.2326678	21.85085477	22.29048838
13.50661475	14.03764956	14.60186131	15.86822818	16.54953392	17.92483736	19.5109009	20.34805237	20.58454167	22.19637462	23.0532186
15.08539834	16.3432109	18.53760668	19.68821266	21.09734591	21.46088109	22.74369911	24.11571756	24.41143433	24.91591262	25.22407808
13.68753534	15.36485386	16.46303757	17.52705112	17.84914043	18.66854213	19.41910279	20.20529032	21.10251303	22.30816275	23.22918935
14.87312764	15.72962855	16.30537959	17.31514318	17.75192705	18.21643446	19.06590519	19.79264017	21.01084206	21.91953723	22.47029613
16.91978721	17.79686763	18.69474301	20.76116567	20.99288194	21.60103775	21.86582688	22.22935373	22.45353612	23.4630839	24.39499851
15.70663563	18.0007559	18.46960568	19.71284641	20.50582527	21.73967091	22.19607942	22.50028049	24.03024962	24.78151088	26.01082011
13.78845438	15.25273768	15.45881165	16.28652792	17.00186011	18.70371152	18.99795568	20.08872059	20.94359519	21.28062804	22.48291439
18.258099	19.89122661	20.69973211	21.28585157	21.94546473	23.25687792	24.43629307	25.45251335	26.19361899	28.82543715	29.49610922
17.48502182	20.17020447	21.98925931	23.47431719	25.01198224	25.48476259	25.73905702	26.25318299	28.03001178	28.98659006	29.32594723
16.4963768	17.49730889	18.98968346	20.034917	20.87026171	21.35518814	23.45581602	24.49822225	24.92167377	25.60058894	25.7674746
15.75606033	16.22008517	18.94879257	19.98268641	21.43947172	23.22800284	23.40625768	23.89369419	24.11038123	24.92194306	25.07252989
17.68272621	18.51089957	19.21089634	19.90813301	21.58315479	22.30212995	22.75009648	24.36597748	26.06461032	26.51947527	26.92251536
15.17895	16.63041328	17.0133685	18.65694291	18.99265988	20.55597182	21.27791097	21.92597898	23.44634256	23.68541471	24.3971786
13.55090502	14.75833388	15.43304986	16.06733741	16.87533985	17.44698278	19.68211096	19.92232368	20.87176888	21.06299848	22.08143017
12.09987271	13.07011494	13.78742954	14.71973487	15.93252705	16.85237435	16.92202902	18.0921888	18.76448862	19.11131359	19.8661034
18.23450576	19.44315577	20.58857353	21.70281273	22.18407485	22.43644514	23.52713505	23.83883205	24.51019453	24.75528924	25.5988776
16.16710918	16.92403543	17.93029195	18.51802878	19.764568	21.87635802	23.40415971	24.7087931	25.91474822	27.08193028	28.0047529
];

J=3;%应力水平数
n=30;%样本个数
l=[11,11,11];%每个应力水平下的监测增量个数

S=[65,85,100];
S0=40;%正常应力水平
SH=100;%容许上限
for k=1:J
    s(k)=(1/(S0+273.15)-1/(S(k)+273.15))/(1/(S0+273.15)-1/(SH+273.15));
end

ff=100;
T1=0:ff:1000;
T2=1000:ff:2000;
T3=2000:ff:3000;

T=0:30:3000;
for j=1:length(T)
    if(T(j)<=1000)
      X_g(j)=u*exp(a*s(1))*T(j)^b;
      Var_g(j)=u^p/lam*exp(a*s(1))*T(j)^b;
    else if (T(j)<=2000)
            tau1=1000;
            v1=0;
            v2=(v1+tau1)*exp(a*(s(1)-s(2))/b);
            X_g(j)=u*exp(a*s(2))*(v2+T(j)-tau1)^b;
            Var_g(j)=u^p/lam*exp(a*s(2))*(v2+T(j)-tau1)^b;
        else
            tau2=2000;
            v3=(v2+tau2-tau1)*exp(a*(s(2)-s(3))/b);
            X_g(j)=u*exp(a*s(3))*(v3+T(j)-tau2)^b;
            Var_g(j)=u^p/lam*exp(a*s(3))*(v3+T(j)-tau2)^b;
        end
    end
end
X_g_L=X_g-norminv(1-0.025,0,1)*sqrt(Var_g);
X_g_U=X_g+norminv(1-0.025,0,1)*sqrt(Var_g);

X1_mean=mean(Y1);
X2_mean=mean(Y2);
X3_mean=mean(Y3);

%%上面是计算数据代码

subplot(1,2,1)
%%下面是绘图代码
%路径图
p1=plot(T1,Y1,'rp:','markersize',3);
hold on
p2=plot(T2,Y2,'bd:','markersize',3);
p3=plot(T3,Y3,'go:','markersize',3);
legend([p1(1) p2(1) p3(1)],{'Deradation paths at 65℃','Deradation paths at 85℃','Deradation paths at 100℃'})
legend('Location','Northwest')

xlabel('Test time (hours)')
ylabel('Light intensity (%)')
axis([0 3000 0 30])
set(gca,'XTick',0:500:3000)
set(gca,'YTick',0:5:30)
legend('boxoff')


subplot(1,2,2)
%估计图
p1=plot(T1,X1_mean,'rp','markersize',5);
hold on
p2=plot(T2,X2_mean,'bd','markersize',5);
p3=plot(T3,X3_mean,'go','markersize',5);
hold on

T11=0:30:1000;
T22=1000:30:2000;
T33=2000:30:3000;
p4=plot(T11,X_g(1:length(T11)),'r-','LineWidth',1.5);
hold on 
p5=plot(T22,X_g((length(T11)):(length(T11)+length(T22))-1),'b-','LineWidth',1.5);
p6=plot(T33,X_g((length(T11)+length(T22))-1:(length(T11)+length(T22)+length(T33))-2),'g-','LineWidth',1.5);
hold on
p7=plot(T11,X_g_L(1:length(T11)),'r:','LineWidth',1.5);
hold on 
p8=plot(T22,X_g_L((length(T11)):(length(T11)+length(T22))-1),'b:','LineWidth',1.5);
p9=plot(T33,X_g_L((length(T11)+length(T22))-1:(length(T11)+length(T22)+length(T33))-2),'g:','LineWidth',1.5);
hold on
p10=plot(T11,X_g_U(1:length(T11)),'r:','LineWidth',1.5);
hold on 
p11=plot(T22,X_g_U((length(T11)):(length(T11)+length(T22))-1),'b:','LineWidth',1.5);
p12=plot(T33,X_g_U((length(T11)+length(T22))-1:(length(T11)+length(T22)+length(T33))-2),'g:','LineWidth',1.5);


legend([p1 p2 p3 p4 p5 p6 p7 p8 p9],{'Sample average at 65℃','Sample average at 85℃','Sample average at 100℃','Estimated mean path at 65℃','Estimated mean path at 85℃','Estimated mean path at 100℃','Estimated 95% CIs at 65℃','Estimated 95% CIs at 85℃','Estimated 95% CIs at 100℃'})
legend('Location','Northwest')

xlabel('Test time (hours)')
ylabel('Light intensity (%)')
axis([0 3000 0 30])
set(gca,'XTick',0:500:3000)
set(gca,'YTick',0:5:30)
legend('boxoff')






