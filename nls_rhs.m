function nls_rhs = nls_rhs(z, ut, dummy, k)

    u = ifft(ut);

    nls_rhs = -(i/2) * (k.^2).*ut + i*fft( (abs(u).^2) .* u );

end