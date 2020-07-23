%% funcion que mide la relacion señal a ruido, [snr,mse]=SNR(Io,IoR);
% inputs:
% Io: imagen original
% IoR: imagen ruidosa
% IoR can contain nans, those points are ignored
% outputs:
% snr= potencia(Io)/potencia(Io-IoR)
% mse= potencia(Io-IoR)

function [snr,mse]=SNR(Io,IoR)

valid_ind = ~isnan(Io(:)) & ~isnan(IoR(:));
if sum(isnan(Io(:)) | isnan(IoR(:)))>0;
    disp('[SNR]> WARNING, Io or IoR contains NaN values')
end

% mean square error
mse = mean((Io(valid_ind)-IoR(valid_ind)).^2);

% we can also use 
% root mean square error: rmse = sqrt(mse)

snr = mean(Io(valid_ind).^2) / mse;