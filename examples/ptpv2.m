function clock_offset = ptp(timestamps)

%%  Clock offset estimation of the IEEE 1588-2008 precision time protocol (PTP)
%   This simple source code implements the PTP offset estimation described in
%   "IEEE Standard for a Precision Clock Synchronization Protocol for 
%	 Networked Measurement and Control Systems,"
%	IEEE 1588-2008 Stadnard, July, 2008.
%	(Working Group: PNCS - Precise Networked Clock Synchronization)
%

%%  Function description:
%   Input: timestamps
%    - timestamps: set of the IEEE 1588 PTP timestamps [t1 C(t2) C(t3) t4 r_s r_d]
%                  where r_s and r_d are the residence time 
%									  of sync and delay-req messages, respectively.
%                  (r_s and r_d were recoded by the CISCO IE-3000-4TC-E)
%
%   Output: clock_offset
%
%%  Implementation
t1 = timestamps(:,1) + timestamps(:,5);		%	timestamps(:,5): residence time of sync message
t2 = timestamps(:,2);
t3 = timestamps(:,3);
t4 = timestamps(:,4) - timestamps(:,6);		%	timestamps(:,6): residence time of delay-req message

packet_delay = 0.5 * ( (t4 - t1) - (t3 - t2) );

if packet_delay <= 0 
	clock_offset = 0;
else
	clock_offset = (t2 - t1) - pd;
end

end
