import * as React from 'react';
import { mount } from 'cypress/react';

function MainApp() {
  return <div>hello world</div>;
}

describe('Health check for Cypress component test', () => {
  it('should mount MainApp', () => {
    mount(<MainApp />);
    cy.get('div').should('include.text', 'hello world');
  });
});
