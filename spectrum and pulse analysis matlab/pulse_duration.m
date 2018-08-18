function T_duration = pulse_duration(T,A,int_point)

bkg = mean(A(1:10));
max_value = max(A-bkg);
processed_data = (A-bkg)/max_value;

[peak,peakloc] = max(processed_data);
% [~,loc1]=min(abs(processed_data(1:peakloc-1)-peak/int_point));
% [~,loc2]=min(abs(processed_data(peakloc:end)-peak/int_point));

[pks,locs] = findpeaks( 1 - (processed_data - peak/int_point).^2 ,'SortStr','descend');

% T_duration = T(peakloc+loc2-1)-T(loc1);
T_duration = abs(T(locs(1)) - T(locs(2)));

% figure(2)
% plot(T,abs(A-peak/2))
