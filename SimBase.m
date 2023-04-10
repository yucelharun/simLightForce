classdef SimBase
    properties
        FieldSize
        Lx, Ly
        Sx, Sy, Fx, Fy
        fps
        FrameNumber
        Dt
        TotalTime 
        Time
        kapha
    end
    methods
        function obj = SimBase(FieldSize, L, fps, FrameNumber, kapha)
            obj.FieldSize = FieldSize;
            obj.Lx = L(2);
            obj.Ly = L(1);
            obj.fps = fps;
            obj.FrameNumber = FrameNumber;
            obj.Sx = -obj.Lx:2*obj.Lx/(obj.FieldSize(2) - 1):obj.Lx;
            obj.Sy = -obj.Ly:2*obj.Ly/(obj.FieldSize(1) - 1):obj.Ly;
            [obj.Fx, obj.Fy] = meshgrid(obj.Sx, obj.Sy); 
            obj.Dt = 1/obj.fps;
            obj.TotalTime = (obj.FrameNumber - 1)*obj.Dt;
            obj.Time = 0:obj.Dt:obj.TotalTime;
            obj.kapha = kapha;
        end
    end
end