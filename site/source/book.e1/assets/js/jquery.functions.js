$(document).ready(function(){

/*
 * 		Top Navigation (Dropdown Menu)
 */


	(function($) {
		$("#nav li.sub").hover(
			function() {
				$(this).find('ul:not(:animated)').fadeIn("fast");
				$('ul:first',this).css('visibility', 'visible');
			},
			function() {
				$(this).find('ul:not(:animated)').fadeOut("fast");
				$('ul:first',this).css('visibility', 'hidden');
			}
		)
	})(jQuery);

           
/*
 * 		Cycle functions (jquery.cycle.js)
 */

	// Add 'scrollVert' functionality for scroll boxe
	(function($) {

		$.fn.cycle.transitions.scrollVert = function($cont, $slides, opts) {
		    $cont.css('overflow','hidden');
		    opts.before.push(function(curr, next, opts, fwd) {
		        $(this).show();
		        var currH = curr.offsetHeight, nextH = next.offsetHeight;
		        opts.cssBefore = fwd ? { top: -nextH } : { top: nextH };
		        opts.animIn.top = 0;
		        opts.animOut.top = fwd ? currH : -currH;
		        $slides.not(curr).css(opts.cssBefore);
		    });
		    opts.cssFirst = { top: 0 };
		    opts.cssAfter = { display: 'none' }
		};

	})(jQuery);


/*
 * 		Latest News Scroller
 */

	$('#scrollNews ul').cycle({ 
	    fx: 'scrollVert',
		speed: 650,
		rev: true,
		timeout: 10000,
		next:   '#scrollNews ol li.next', 
	    prev:   '#scrollNews ol li.previous'
	});


/*
 * 		Spotlight Scroller
 */

	$('#scrollSpotlight ul').cycle({ 
	    fx: 'scrollVert',
		speed: 850,
		rev: true,
		timeout: 10000,
		next:   '#scrollSpotlight ol li.next', 
	    prev:   '#scrollSpotlight ol li.previous'
	});


/*
 * 		Contact Form Validation
 */


	$('#contactform').submit(function() {
	
		// Disable the submit button
		$('#contactform input[type=submit]')
			.attr('value', 'Sending message…')
			.attr('disabled', 'disabled');
	
		// AJAX POST request
		$.post(
			$(this).attr('action'),
			{
				name:$('#name').val(),
				email:$('#email').val(),
				message:$('#message').val()
			},
			function(errors) {
				// No errors
				if (errors == null) {
					$('#contactform')
						.hide()
						.html('<h3>Thank you</h3><p>Your message has been sent.</p>')
						.show();
				}
	
				// Errors
				else {
					// Re-enable the submit button
					$('#contactform input[type=submit]')
						.removeAttr('disabled')
						.attr('value', 'Send your Question');
	
					// Technical server problem, the email could not be sent
					if (errors.server != null) {
						alert(errors.server);
						return false;
					}
	
					// Empty the errorbox and reset the error alerts
					$('#contactform .errorbox').html('<ul></ul>').show();
					$('#contactform li').removeClass('alert');
	
					// Loop over the errors, mark the corresponding input fields,
					// and add the error messages to the errorbox.
					for (field in errors) {
						if (errors[field] != null) {
							$('#' + field).parent('li').addClass('alert');
							$('#contactform .errorbox ul').append('<li>' + errors[field] + '</li>');
						}
					}
				}
			},
			'json'
		);
	
		// Prevent non-AJAX form submission
		return false;
	});

});
