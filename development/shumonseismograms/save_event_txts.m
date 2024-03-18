clear all; close all; clc;
load('Events_VR_DATA.mat');
num_of_events = length(Events_sbi_time);
for event_i = 1 : num_of_events
    %acess the data
    event_i_time = Events_sbi_time{event_i};
    event_i_vxvy_bot = Events_sbivel_b{event_i};
    %extract vx
    event_i_vx_bot = event_i_vxvy_bot(1:2:end,:);
    %extract vy
    event_i_vy_bot = event_i_vxvy_bot(2:2:end,:);
    %save data
    writematrix(event_i_time,"./events_data/time_event"+string(event_i)+".txt")
    writematrix(event_i_vx_bot,"./events_data/vx_event"+string(event_i)+".txt")
    writematrix(event_i_vy_bot,"./events_data/vy_event"+string(event_i)+".txt")
end
