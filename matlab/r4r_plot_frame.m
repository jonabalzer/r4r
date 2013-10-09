function [ output_args ] = r4r_plot_frame( frame, scale, handle, annotation )

    figure(handle);
    
    hold on;

    text(frame(1,4),frame(2,4),frame(3,4),annotation);
    
    v = [1;0;0;0]; rv = frame * v; rvs = rv(1:3) * scale;
    quiver3(frame(1,4),frame(2,4),frame(3,4),rvs(1),rvs(2),rvs(3),0,'r');
    v = [0;1;0;0]; rv = frame * v; rvs = rv(1:3) * scale;
    quiver3(frame(1,4),frame(2,4),frame(3,4),rvs(1),rvs(2),rvs(3),0,'g');
    v = [0;0;1;0]; rv = frame * v; rvs = rv(1:3) * scale;
    quiver3(frame(1,4),frame(2,4),frame(3,4),rvs(1),rvs(2),rvs(3),0,'b');
    
    hold off;

    axis equal;

end

