 %% face detection
clear
        VidFrame = imread('im.png');
        % Create a cascade detector object.
        faceDetector = vision.CascadeObjectDetector();

        bbox         = step(faceDetector, VidFrame);
       
        if ~isempty(bbox)
             % Draw the returned bounding box around the detected face.
            VidFrame = insertShape(VidFrame, 'Rectangle', bbox);

            VidFrame = imcrop(VidFrame, bbox(1, :));
        end
      
        %position for optional face detection/tracking - originally specified in reference as a manual segmentation.
       imshow(VidFrame);
        
        