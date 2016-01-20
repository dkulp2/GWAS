<?php

// CONFIGURATION --------------------------------------------------------------

// This is the email where the contact mails will be sent to.
$config['recipient'] = 'foulkes@schoolph.umass.edu';

// This is the subject line for contact emails.
// The variable %name% will be replaced with the name of the sender.
$config['subject'] = 'Contact message from %name%';

// These are the messages displayed in case of form errors.
$config['errors'] = array
(
	'no_name'       => 'Please enter your name',
	'no_email'      => 'Your Email adresse is required.',
	'invalid_email' => 'You entered an invalid email address.',
	'no_message'    => 'Please, include your message.',
);

// END OF CONFIGURATION -------------------------------------------------------


// Ignore non-POST requests
if ( ! $_POST)
	exit('Nothing to see here. Please go back to the site.');

// Was this an AJAX request or not?
$ajax = (isset($_SERVER['HTTP_X_REQUESTED_WITH']) AND strtolower($_SERVER['HTTP_X_REQUESTED_WITH']) == 'xmlhttprequest');

// Set the correct HTTP headers
header('Content-Type: text/'.($ajax ? 'plain' : 'html').'; charset=utf-8');

// Extract and trim contactform values
$name    = isset($_POST['name']) ? trim($_POST['name']) : '';
$email   = isset($_POST['email']) ? trim($_POST['email']) : '';
$message = isset($_POST['message']) ? trim($_POST['message']) : '';

// Take care of magic quotes if needed (you really should have them disabled)
set_magic_quotes_runtime(0);
if (get_magic_quotes_gpc())
{
	$name    = stripslashes($name);
	$email   = stripslashes($email);
	$message = stripslashes($message);
}

// Initialize the errors array which will also be sent back as a JSON object
$errors = NULL;

// Validate name
if ($name == '' || strpos($name, "\r") || strpos($name, "\n"))
{
	$errors['name'] = $config['errors']['no_name'];
}

// Validate email
if ($email == '')
{
	$errors['email'] = $config['errors']['no_email'];
}
elseif ( ! preg_match('/^[-_a-z0-9\'+*$^&%=~!?{}]++(?:\.[-_a-z0-9\'+*$^&%=~!?{}]+)*+@(?:(?![-.])[-a-z0-9.]+(?<![-.])\.[a-z]{2,6}|\d{1,3}(?:\.\d{1,3}){3})(?::\d++)?$/iD', $email))
{
	$errors['email'] = $config['errors']['invalid_email'];
}

// Validate message
if ($message == '')
{
	$errors['message'] = $config['errors']['no_message'];
}

// Validation succeeded
if (empty($errors))
{
	// Prepare subject line
	$subject = str_replace('%name%', $name, $config['subject']);

	// Additional mail headers
	$headers  = 'Content-Type: text/plain; charset=utf-8'."\r\n";
	$headers .= 'From: '.$email;

	// Send the mail
	if ( ! mail($config['recipient'], $subject, $message, $headers))
	{
		$errors['server'] = 'There seems to be a technical problem with our server. We are sorry. '.
		                    'Could you mail your message directly at '.$config['recipient'].'? Thank you.';
	}
}

if ($ajax)
{
	// Output the possible errors as a JSON object
	echo json_encode($errors);
}
else
{
	// Show a simple HTML feedback message in case of non-javascript support
	if (empty($errors))
	{
		echo '<h1>Thank you</h1>';
		echo '<p>Your message has been sent.</p>';
	}
	else
	{
		echo '<h1>Oops!</h1>';
		echo '<p>Please go back and fix the following errors:</p>';
		echo '<ul><li>';
		echo implode('</li><li>', $errors);
		echo '</li></ul>';
	}
}