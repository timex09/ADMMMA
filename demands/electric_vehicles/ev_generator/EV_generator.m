function [EV_daily_occurence,EV_daily_demand,EV_uncontrolled_charging_profile,SOC_uncontrolled_charging_profile]=EV_generator(number_profiles,charging_power, year, charging_eff,battery_capacity,e_demand_specific, arrival_prob,departure_prob,daily_distance_prob,distance_prob_factor,velocity_prob)

anz = 365;

startDate =datenum(['Mon 01.01.' num2str(year) ],'ddd dd.mm.yyyy');

d = zeros(anz, 1);

d(1)= startDate;
for k = 1 : anz
   d(k+1)= startDate + k;
end

datum=datestr(d,'ddd dd.mm.yyyy');

EV_daily_occurence=zeros(anz*size(departure_prob,1),number_profiles);
EV_daily_demand=zeros(anz,number_profiles);
EV_uncontrolled_charging_profile=zeros(anz*size(departure_prob,1),number_profiles);
SOC_uncontrolled_charging_profile=zeros(anz*size(departure_prob,1),number_profiles);

%%
parfor jj=1:number_profiles
%disp(jj);
    occurence_result=ones(96,anz+1);
    E_daily_demand_result=zeros(anz,1);
    EV_uncontrolled_charging_result=zeros(96,anz+1);
    Soc_Bat_uncontrolled_result=zeros(96,anz+1);
    battery_capacity_result=battery_capacity;

    for ii=1:anz

        if strcmp(datum(ii,1:3),'Sat') 

            [arrival]=randi_probability(arrival_prob(:,2),4*1,1:96);
            y=accumarray(arrival(:),1);
            [~,I]=max(y);
            arrival=I;
            [departure]=randi_probability(departure_prob(:,2),4*1,1:96);
            y=accumarray(departure(:),1);
            [~,I]=max(y);
            departure=I;
            distance=mean(daily_distance_prob(:,6))+std(daily_distance_prob(:,6)).*randn(1,1);
            a=find(~xor(distance_prob_factor(1,:)>distance, distance_prob_factor(2,:)<distance))-1;
            b=find(~xor(distance_prob_factor(1,:)>distance, distance_prob_factor(2,:)<distance))+1;
                if a==0 
                    a=1;
                    b=2;
                end
                
                if b==7
                    a=5;
                    b=6;
                end                
            distance=mean(1/10*round(10*((distance_prob_factor(2,a:b)-distance_prob_factor(1,a:b))*rand(1,1)+distance_prob_factor(1,a:b))));[velocity]=randi_probability(velocity_prob,2,1:120);
            [velocity]=randi_probability(velocity_prob,2,1:120);
            velocity=mean(velocity);

            while arrival<departure || abs(arrival-departure)<ceil(distance/velocity/0.25)

                [arrival]=randi_probability(arrival_prob(:,2),4*1,1:96);
                y=accumarray(arrival(:),1);
                [~,I]=max(y);
                arrival=I;
                [departure]=randi_probability(departure_prob(:,2),4*1,1:96);
                y=accumarray(departure(:),1);
                [~,I]=max(y);
                departure=I;

            end

            occurence_result(departure:arrival-1,ii)=zeros(arrival-departure,1);
            E_daily_demand_result(ii,1)=distance*e_demand_specific;

            load_profile_=[charging_power*ones(floor(E_daily_demand_result(ii,1)*4/(charging_power*charging_eff)),1); mod(E_daily_demand_result(ii,1)*4/(charging_power*charging_eff),1)*charging_power];
            load_profile=zeros(2*96,1);
            load_profile(arrival:arrival+size(load_profile_,1)-1)=load_profile_;
            load_profile(~[occurence_result(:,ii); occurence_result(:,ii+1)])=0;
            EV_uncontrolled_charging_result(:,ii:ii+1)=reshape(load_profile,96,2);
            
            soc_profile=charging_eff*load_profile/4;
            soc_profile(arrival-1)=battery_capacity_result-E_daily_demand_result(ii);
            soc_profile(arrival-1:arrival+size(load_profile_,1)-1)=cumsum(soc_profile((arrival-1:arrival+size(load_profile_,1)-1)));
            battery_capacity_result=soc_profile(arrival+size(load_profile_,1)-1);
            Soc_Bat_uncontrolled_result(:,ii:ii+1)=reshape(soc_profile,96,2);

        elseif strcmp(datum(ii,1:3),'Sun') 

            [arrival]=randi_probability(arrival_prob(:,3),4*1,1:96);
            y=accumarray(arrival(:),1);
            [~,I]=max(y);
            arrival=I;
            [departure]=randi_probability(departure_prob(:,3),4*1,1:96);
            y=accumarray(departure(:),1);
            [~,I]=max(y);
            departure=I;
            distance=mean(daily_distance_prob(:,7))+std(daily_distance_prob(:,7)).*randn(1,1);
            a=find(~xor(distance_prob_factor(1,:)>distance, distance_prob_factor(2,:)<distance))-1;
            b=find(~xor(distance_prob_factor(1,:)>distance, distance_prob_factor(2,:)<distance))+1;
                if a==0 
                    a=1;
                    b=2;
                end
                
                if b==7
                    a=5;
                    b=6;
                end                
            distance=mean(1/10*round(10*((distance_prob_factor(2,a:b)-distance_prob_factor(1,a:b))*rand(1,1)+distance_prob_factor(1,a:b))));[velocity]=randi_probability(velocity_prob,2,1:120);
            velocity=mean(velocity);

            while arrival<departure || abs(arrival-departure)<ceil(distance/velocity/0.25)

                [arrival]=randi_probability(arrival_prob(:,3),4*1,1:96);
                y=accumarray(arrival(:),1);
                [~,I]=max(y);
                arrival=I;
                [departure]=randi_probability(departure_prob(:,3),4*1,1:96);
                y=accumarray(departure(:),1);
                [~,I]=max(y);
                departure=I;

            end

            occurence_result(departure:arrival-1,ii)=zeros(arrival-departure,1);
            E_daily_demand_result(ii,1)=distance*e_demand_specific;

            load_profile_=[charging_power*ones(floor(E_daily_demand_result(ii,1)*4/(charging_power*charging_eff)),1); mod(E_daily_demand_result(ii,1)*4/(charging_power*charging_eff),1)*charging_power];
            load_profile=zeros(2*96,1);
            load_profile(arrival:arrival+size(load_profile_,1)-1)=load_profile_;
            load_profile(~[occurence_result(:,ii); occurence_result(:,ii+1)])=0;
            EV_uncontrolled_charging_result(:,ii:ii+1)=reshape(load_profile,96,2);
            
            soc_profile=charging_eff*load_profile/4;
            soc_profile(arrival-1)=battery_capacity_result-E_daily_demand_result(ii);
            soc_profile(arrival-1:arrival+size(load_profile_,1)-1)=cumsum(soc_profile((arrival-1:arrival+size(load_profile_,1)-1)));
            battery_capacity_result=soc_profile(arrival+size(load_profile_,1)-1);
            Soc_Bat_uncontrolled_result(:,ii:ii+1)=reshape(soc_profile,96,2);
                                   

        elseif strcmp(datum(ii,1:3),'Mon') 

            [arrival]=randi_probability(arrival_prob(:,1),4*5,1:96);
            y=accumarray(arrival(:),1);
            [~,I]=max(y);
            arrival=I;
            [departure]=randi_probability(departure_prob(:,1),4*5,1:96);
            y=accumarray(departure(:),1);
            [~,I]=max(y);
            departure=I;
            distance=mean(daily_distance_prob(:,1))+std(daily_distance_prob(:,1)).*randn(1,1);
            a=find(~xor(distance_prob_factor(1,:)>distance, distance_prob_factor(2,:)<distance))-1;
            b=find(~xor(distance_prob_factor(1,:)>distance, distance_prob_factor(2,:)<distance))+1;
                if a==0 
                    a=1;
                    b=2;
                end
                
                if b==7
                    a=5;
                    b=6;
                end                
            distance=mean(1/10*round(10*((distance_prob_factor(2,a:b)-distance_prob_factor(1,a:b))*rand(1,1)+distance_prob_factor(1,a:b))));
            [velocity]=randi_probability(velocity_prob,2,1:120);
            velocity=mean(velocity);

            while arrival<departure || abs(arrival-departure)<ceil(distance/velocity/0.25)

                [arrival]=randi_probability(arrival_prob(:,1),4*5,1:96);
                y=accumarray(arrival(:),1);
                [~,I]=max(y);
                arrival=I;
                [departure]=randi_probability(departure_prob(:,1),4*5,1:96);
                y=accumarray(departure(:),1);
                [~,I]=max(y);
                departure=I;

            end

            occurence_result(departure:arrival-1,ii)=zeros(arrival-departure,1);
            E_daily_demand_result(ii,1)=distance*e_demand_specific;

            load_profile_=[charging_power*ones(floor(E_daily_demand_result(ii,1)*4/(charging_power*charging_eff)),1); mod(E_daily_demand_result(ii,1)*4/(charging_power*charging_eff),1)*charging_power];
            load_profile=zeros(2*96,1);
            load_profile(arrival:arrival+size(load_profile_,1)-1)=load_profile_;
            load_profile(~[occurence_result(:,ii); occurence_result(:,ii+1)])=0;
            EV_uncontrolled_charging_result(:,ii:ii+1)=reshape(load_profile,96,2);
            
            soc_profile=charging_eff*load_profile/4;
            soc_profile(arrival-1)=battery_capacity_result-E_daily_demand_result(ii);
            soc_profile(arrival-1:arrival+size(load_profile_,1)-1)=cumsum(soc_profile((arrival-1:arrival+size(load_profile_,1)-1)));
            battery_capacity_result=soc_profile(arrival+size(load_profile_,1)-1);
            Soc_Bat_uncontrolled_result(:,ii:ii+1)=reshape(soc_profile,96,2);

        elseif strcmp(datum(ii,1:3),'Tue') 

            [arrival]=randi_probability(arrival_prob(:,1),4*5,1:96);
            y=accumarray(arrival(:),1);
            [~,I]=max(y);
            arrival=I;
            [departure]=randi_probability(departure_prob(:,1),4*5,1:96);
            y=accumarray(departure(:),1);
            [~,I]=max(y);
            departure=I;
            distance=mean(daily_distance_prob(:,2))+std(daily_distance_prob(:,2)).*randn(1,1);
            a=find(~xor(distance_prob_factor(1,:)>distance, distance_prob_factor(2,:)<distance))-1;
            b=find(~xor(distance_prob_factor(1,:)>distance, distance_prob_factor(2,:)<distance))+1;
                if a==0 
                    a=1;
                    b=2;
                end
                
                if b==7
                    a=5;
                    b=6;
                end                
            distance=mean(1/10*round(10*((distance_prob_factor(2,a:b)-distance_prob_factor(1,a:b))*rand(1,1)+distance_prob_factor(1,a:b))));[velocity]=randi_probability(velocity_prob,2,1:120);
            [velocity]=randi_probability(velocity_prob,2,1:120);
            velocity=mean(velocity);

            while arrival<departure || abs(arrival-departure)<ceil(distance/velocity/0.25)

                [arrival]=randi_probability(arrival_prob(:,1),4*5,1:96);
                y=accumarray(arrival(:),1);
                [~,I]=max(y);
                arrival=I;
                [departure]=randi_probability(departure_prob(:,1),4*5,1:96);
                y=accumarray(departure(:),1);
                [~,I]=max(y);
                departure=I;

            end 

            occurence_result(departure:arrival-1,ii)=zeros(arrival-departure,1);
            E_daily_demand_result(ii,1)=distance*e_demand_specific;

            load_profile_=[charging_power*ones(floor(E_daily_demand_result(ii,1)*4/(charging_power*charging_eff)),1); mod(E_daily_demand_result(ii,1)*4/(charging_power*charging_eff),1)*charging_power];
            load_profile=zeros(2*96,1);
            load_profile(arrival:arrival+size(load_profile_,1)-1)=load_profile_;
            load_profile(~[occurence_result(:,ii); occurence_result(:,ii+1)])=0;
            EV_uncontrolled_charging_result(:,ii:ii+1)=reshape(load_profile,96,2);
            
            soc_profile=charging_eff*load_profile/4;
            soc_profile(arrival-1)=battery_capacity_result-E_daily_demand_result(ii);
            soc_profile(arrival-1:arrival+size(load_profile_,1)-1)=cumsum(soc_profile((arrival-1:arrival+size(load_profile_,1)-1)));
            battery_capacity_result=soc_profile(arrival+size(load_profile_,1)-1);
            Soc_Bat_uncontrolled_result(:,ii:ii+1)=reshape(soc_profile,96,2);

        elseif strcmp(datum(ii,1:3),'Wed') 

            [arrival]=randi_probability(arrival_prob(:,1),4*5,1:96);
            y=accumarray(arrival(:),1);
            [~,I]=max(y);
            arrival=I;
            [departure]=randi_probability(departure_prob(:,1),4*5,1:96);
            y=accumarray(departure(:),1);
            [~,I]=max(y);
            departure=I;
            distance=mean(daily_distance_prob(:,3))+std(daily_distance_prob(:,3)).*randn(1,1);
            a=find(~xor(distance_prob_factor(1,:)>distance, distance_prob_factor(2,:)<distance))-1;
            b=find(~xor(distance_prob_factor(1,:)>distance, distance_prob_factor(2,:)<distance))+1;
                if a==0 
                    a=1;
                    b=2;
                end
                
                if b==7
                    a=5;
                    b=6;
                end                
            distance=mean(1/10*round(10*((distance_prob_factor(2,a:b)-distance_prob_factor(1,a:b))*rand(1,1)+distance_prob_factor(1,a:b))));[velocity]=randi_probability(velocity_prob,2,1:120);
            [velocity]=randi_probability(velocity_prob,2,1:120);
            velocity=mean(velocity);

            while arrival<departure || abs(arrival-departure)<ceil(distance/velocity/0.25)

                [arrival]=randi_probability(arrival_prob(:,1),4*5,1:96);
                y=accumarray(arrival(:),1);
                [~,I]=max(y);
                arrival=I;
                [departure]=randi_probability(departure_prob(:,1),4*5,1:96);
                y=accumarray(departure(:),1);
                [~,I]=max(y);
                departure=I;

            end   

            occurence_result(departure:arrival-1,ii)=zeros(arrival-departure,1);
            E_daily_demand_result(ii,1)=distance*e_demand_specific;

            load_profile_=[charging_power*ones(floor(E_daily_demand_result(ii,1)*4/(charging_power*charging_eff)),1); mod(E_daily_demand_result(ii,1)*4/(charging_power*charging_eff),1)*charging_power];
            load_profile=zeros(2*96,1);
            load_profile(arrival:arrival+size(load_profile_,1)-1)=load_profile_;
            load_profile(~[occurence_result(:,ii); occurence_result(:,ii+1)])=0;
            EV_uncontrolled_charging_result(:,ii:ii+1)=reshape(load_profile,96,2);
            
            soc_profile=charging_eff*load_profile/4;
            soc_profile(arrival-1)=battery_capacity_result-E_daily_demand_result(ii);
            soc_profile(arrival-1:arrival+size(load_profile_,1)-1)=cumsum(soc_profile((arrival-1:arrival+size(load_profile_,1)-1)));
            battery_capacity_result=soc_profile(arrival+size(load_profile_,1)-1);
            Soc_Bat_uncontrolled_result(:,ii:ii+1)=reshape(soc_profile,96,2);

        elseif strcmp(datum(ii,1:3),'Thu') 

            [arrival]=randi_probability(arrival_prob(:,1),4*5,1:96);
            y=accumarray(arrival(:),1);
            [~,I]=max(y);
            arrival=I;
            [departure]=randi_probability(departure_prob(:,1),4*5,1:96);
            y=accumarray(departure(:),1);
            [~,I]=max(y);
            departure=I;
            distance=mean(daily_distance_prob(:,4))+std(daily_distance_prob(:,4)).*randn(1,1);
            a=find(~xor(distance_prob_factor(1,:)>distance, distance_prob_factor(2,:)<distance))-1;
            b=find(~xor(distance_prob_factor(1,:)>distance, distance_prob_factor(2,:)<distance))+1;
                if a==0 
                    a=1;
                    b=2;
                end
                
                if b==7
                    a=5;
                    b=6;
                end                
            distance=mean(1/10*round(10*((distance_prob_factor(2,a:b)-distance_prob_factor(1,a:b))*rand(1,1)+distance_prob_factor(1,a:b))));[velocity]=randi_probability(velocity_prob,2,1:120);
            [velocity]=randi_probability(velocity_prob,2,1:120);
            velocity=mean(velocity);

            while arrival<departure || abs(arrival-departure)<ceil(distance/velocity/0.25)

                [arrival]=randi_probability(arrival_prob(:,1),4*5,1:96);
                y=accumarray(arrival(:),1);
                [~,I]=max(y);
                arrival=I;
                [departure]=randi_probability(departure_prob(:,1),4*5,1:96);
                y=accumarray(departure(:),1);
                [~,I]=max(y);
                departure=I;

            end   

            occurence_result(departure:arrival-1,ii)=zeros(arrival-departure,1);
            E_daily_demand_result(ii,1)=distance*e_demand_specific;

            load_profile_=[charging_power*ones(floor(E_daily_demand_result(ii,1)*4/(charging_power*charging_eff)),1); mod(E_daily_demand_result(ii,1)*4/(charging_power*charging_eff),1)*charging_power];
            load_profile=zeros(2*96,1);
            load_profile(arrival:arrival+size(load_profile_,1)-1)=load_profile_;
            load_profile(~[occurence_result(:,ii); occurence_result(:,ii+1)])=0;
            EV_uncontrolled_charging_result(:,ii:ii+1)=reshape(load_profile,96,2);
            
            soc_profile=charging_eff*load_profile/4;
            soc_profile(arrival-1)=battery_capacity_result-E_daily_demand_result(ii);
            soc_profile(arrival-1:arrival+size(load_profile_,1)-1)=cumsum(soc_profile((arrival-1:arrival+size(load_profile_,1)-1)));
            battery_capacity_result=soc_profile(arrival+size(load_profile_,1)-1);
            Soc_Bat_uncontrolled_result(:,ii:ii+1)=reshape(soc_profile,96,2);

        elseif strcmp(datum(ii,1:3),'Fri') 

            [arrival]=randi_probability(arrival_prob(:,1),4*5,1:96);
            y=accumarray(arrival(:),1);
            [~,I]=max(y);
            arrival=I;
            [departure]=randi_probability(departure_prob(:,1),4*5,1:96);
            y=accumarray(departure(:),1);
            [~,I]=max(y);
            departure=I;
            distance=mean(daily_distance_prob(:,5))+std(daily_distance_prob(:,5)).*randn(1,1);
            a=find(~xor(distance_prob_factor(1,:)>distance, distance_prob_factor(2,:)<distance))-1;
            b=find(~xor(distance_prob_factor(1,:)>distance, distance_prob_factor(2,:)<distance))+1;
                if a==0 
                    a=1;
                    b=2;
                end
                
                if b==7
                    a=5;
                    b=6;
                end                
            distance=mean(1/10*round(10*((distance_prob_factor(2,a:b)-distance_prob_factor(1,a:b))*rand(1,1)+distance_prob_factor(1,a:b))));[velocity]=randi_probability(velocity_prob,2,1:120);
            [velocity]=randi_probability(velocity_prob,2,1:120);
            velocity=mean(velocity);

            while arrival<departure || abs(arrival-departure)<ceil(distance/velocity/0.25)

                [arrival]=randi_probability(arrival_prob(:,1),4*5,1:96);
                y=accumarray(arrival(:),1);
                [~,I]=max(y);
                arrival=I;
                [departure]=randi_probability(departure_prob(:,1),4*5,1:96);
                y=accumarray(departure(:),1);
                [~,I]=max(y);
                departure=I;

            end   

            occurence_result(departure:arrival-1,ii)=zeros(arrival-departure,1);
            E_daily_demand_result(ii,1)=distance*e_demand_specific;

            load_profile_=[charging_power*ones(floor(E_daily_demand_result(ii,1)*4/(charging_power*charging_eff)),1); mod(E_daily_demand_result(ii,1)*4/(charging_power*charging_eff),1)*charging_power];
            load_profile=zeros(2*96,1);
            load_profile(arrival:arrival+size(load_profile_,1)-1)=load_profile_;
            load_profile(~[occurence_result(:,ii); occurence_result(:,ii+1)])=0;
            EV_uncontrolled_charging_result(:,ii:ii+1)=reshape(load_profile,96,2);
            
            soc_profile=charging_eff*load_profile/4;
            soc_profile(arrival-1)=battery_capacity_result-E_daily_demand_result(ii);
            soc_profile(arrival-1:arrival+size(load_profile_,1)-1)=cumsum(soc_profile((arrival-1:arrival+size(load_profile_,1)-1)));
            battery_capacity_result=soc_profile(arrival+size(load_profile_,1)-1);
            Soc_Bat_uncontrolled_result(:,ii:ii+1)=reshape(soc_profile,96,2);
            
        end      

    end   

occurence_result(:,366)=[];
occurence_result=occurence_result(:);
EV_uncontrolled_charging_result(:,366)=[];
EV_uncontrolled_charging_result=EV_uncontrolled_charging_result(:);
Soc_Bat_uncontrolled_result(:,366)=[];
Soc_Bat_uncontrolled_result=Soc_Bat_uncontrolled_result(:);
    
    
EV_daily_occurence(:,jj)=occurence_result;
EV_daily_demand(:,jj)=E_daily_demand_result;
EV_uncontrolled_charging_profile(:,jj)=EV_uncontrolled_charging_result;
SOC_uncontrolled_charging_profile(:,jj)=Soc_Bat_uncontrolled_result;

end

end